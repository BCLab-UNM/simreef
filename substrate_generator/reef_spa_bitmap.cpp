// reef_spa_substrate.cpp
// Power-law patch sizes + Perlin/fBm blobs with warp + tunable params

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>

static const int DX4[4] = {1,-1,0,0};
static const int DY4[4] = {0,0,1,-1};
inline bool inside(int x,int y,int W,int H){return x>=0&&x<W&&y>=0&&y<H;}

// ---------------------
// Connected component cleanup
// ---------------------
struct Component {int id; int color; std::vector<std::pair<int,int>> pixels;};

void enforce_min_patch(std::vector<std::vector<int>>& img,int W,int H,int min_patch){
    if(min_patch<=0) return;
    bool changed=true;
    while(changed){
        changed=false;
        std::vector<std::vector<int>> labels(H, std::vector<int>(W,-1));
        int next_id=0;
        std::unordered_map<int,Component> comps;

        for(int y=0;y<H;y++){
            for(int x=0;x<W;x++){
                if(labels[y][x]!=-1) continue;
                int col=img[y][x];
                std::queue<std::pair<int,int>> q;
                q.push({x,y});
                labels[y][x]=next_id;
                comps[next_id]={next_id,col,{}};
                while(!q.empty()){
                    auto [cx,cy]=q.front(); q.pop();
                    comps[next_id].pixels.push_back({cx,cy});
                    for(int k=0;k<4;k++){
                        int nx=cx+DX4[k], ny=cy+DY4[k];
                        if(inside(nx,ny,W,H) && labels[ny][nx]==-1 && img[ny][nx]==col){
                            labels[ny][nx]=next_id;
                            q.push({nx,ny});
                        }
                    }
                }
                next_id++;
            }
        }

        for(auto &kv : comps){
            auto &comp = kv.second;
            if((int)comp.pixels.size() < min_patch){
                std::unordered_map<int,int> border;
                for(auto &pt : comp.pixels){
                    int x=pt.first, y=pt.second;
                    for(int k=0;k<4;k++){
                        int nx=x+DX4[k], ny=y+DY4[k];
                        if(inside(nx,ny,W,H) && labels[ny][nx]!=comp.id){
                            border[ img[ny][nx] ]++;
                        }
                    }
                }
                int best=comp.color, bc=-1;
                for(auto &b : border){ if(b.second>bc){ bc=b.second; best=b.first; } }
                for(auto &pt : comp.pixels) img[pt.second][pt.first]=best;
                changed=true;
            }
        }
    }
}

// ---------------------
// BMP writer
// ---------------------
#pragma pack(push,1)
struct BMPHeader {
    uint16_t bfType{0x4D42};
    uint32_t bfSize{0};
    uint16_t bfReserved1{0};
    uint16_t bfReserved2{0};
    uint32_t bfOffBits{54};
    uint32_t biSize{40};
    int32_t  biWidth{0};
    int32_t  biHeight{0};
    uint16_t biPlanes{1};
    uint16_t biBitCount{24};
    uint32_t biCompression{0};
    uint32_t biSizeImage{0};
    int32_t  biXPelsPerMeter{2835};
    int32_t  biYPelsPerMeter{2835};
    uint32_t biClrUsed{0};
    uint32_t biClrImportant{0};
};
#pragma pack(pop)

void writeBMP(const std::string& fn,const std::vector<std::vector<int>>& img,int W,int H){
    BMPHeader hdr;
    hdr.biWidth = W; hdr.biHeight = H;
    int rowSize = ((W*3+3)&~3);
    hdr.biSizeImage = rowSize * H;
    hdr.bfSize = hdr.bfOffBits + hdr.biSizeImage;

    std::ofstream out(fn, std::ios::binary);
    out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

    
    //uint8_t pal[4][3] = {
    //{  0,   0, 255}, // Red
    //{  0, 255,   0}, // Green
    //{255,   0,   0}, // Blue
    //{  0, 255, 255}  // Yellow
    //};
    

    // RGB
    //uint8_t pal[4][3] = {
    //			 {255, 255,   0}, // 0 - Sand (Yellow)
    //			 {  0, 255,   0}, // 1 - Algae on Sand (Green)
    //			 {255,   0,   0}, // 2 - Coral (Red)
    //			 {139,  69,  19}  // 3 - Algae on Coral (Brown)
    //};
    
    
    // BGR
    //uint8_t pal[4][3] = {
    //			 {  0, 255, 255}, // 0 - Sand (Yellow)   (B=0, G=255, R=255)
    //			 {  0, 255,   0}, // 1 - Algae on Sand (Green)
    //			 {  0,   0, 255}, // 2 - Coral (Red)
    //			 { 19,  69, 139}  // 3 - Algae on Coral (Brown)
    //};

    // BGR palette (with olive tones for algae on coral)
    uint8_t pal[4][3] = {
			 {0, 255, 255}, // 0 - Sand (light yellow)
			 {  0, 255,   0}, // 1 - Algae on Sand (vibrant green)
			 {0, 0,   255}, // 2 - Coral with no algae
			 { 0, 180,  165}  // 3 - Algae on Coral (dark olive green)
    };

    
    std::vector<uint8_t> row(rowSize,0);
    for(int y=0;y<H;y++){
        int yy = H-1-y;
        for(int x=0;x<W;x++){
            int c = img[yy][x] & 3;
            row[x*3+0] = pal[c][0];
            row[x*3+1] = pal[c][1];
            row[x*3+2] = pal[c][2];
        }
        out.write(reinterpret_cast<char*>(row.data()), rowSize);
    }
}

// ---------------------
// Perlin noise + fBm
// ---------------------
double fade(double t){ return t*t*t*(t*(t*6-15)+10); }
double lerp(double a,double b,double t){ return a+(b-a)*t; }

struct Perlin {
    int p[512];
    Perlin(unsigned seed=0){
        std::iota(p,p+256,0);
        std::mt19937 rng(seed);
        std::shuffle(p,p+256,rng);
        for(int i=0;i<256;i++) p[256+i]=p[i];
    }
    double grad(int hash,double x,double y){
        int h=hash&7;
        double u=h<4?x:y;
        double v=h<4?y:x;
        return ((h&1)?-u:u)+((h&2)?-v:v);
    }
    double noise(double x,double y) {
        int X=(int)floor(x)&255;
        int Y=(int)floor(y)&255;
        x-=floor(x); y-=floor(y);
        double u=fade(x), v=fade(y);
        int A=p[X]+Y, B=p[X+1]+Y;
        return lerp( lerp(grad(p[A],x,y), grad(p[B],x-1,y), u),
                     lerp(grad(p[A+1],x,y-1), grad(p[B+1],x-1,y-1), u), v);
    }
};

// fBm with domain warp
double fBmWarp(Perlin& perlin,double x,double y,
               double base_freq,int octaves,double persistence,
               double lacunarity,double warp_strength){
    double wx = x + warp_strength * perlin.noise(x*0.001, y*0.001);
    double wy = y + warp_strength * perlin.noise(x*0.001+100, y*0.001+100);

    double val=0, amp=1, freq=base_freq;
    for(int o=0;o<octaves;o++){
        val += amp * perlin.noise(wx*freq, wy*freq);
        freq *= lacunarity;
        amp *= persistence;
    }
    return val;
}

// ---------------------
// Power-law sampling
// ---------------------
double alpha_to_tau(double alpha){
    double a = std::max(0.2, alpha);
    return 1.0 + 1.0 / a;
}
int sample_pareto_size(std::mt19937& rng, int min_patch, double tau, int max_cap){
    if(min_patch<=0) min_patch=1;
    std::uniform_real_distribution<double> U(0.0,1.0);
    double u = std::max(1e-12,std::min(1-1e-12,U(rng)));
    double x = (double)min_patch / std::pow(1.0-u,1.0/tau);
    long long s=(long long)std::llround(x);
    if(s<min_patch) s=min_patch;
    if(max_cap>0 && s>max_cap) s=max_cap;
    return (int)s;
}

// ---------------------
// Blob growth with compactness
// ---------------------
struct Node {int x,y; double d;};
struct NodeCmp { bool operator()(const Node&a,const Node&b) const {return a.d>b.d;} };

void grow_blob_perlin(std::vector<std::vector<int>>& img,
                      int W,int H,
                      int seedx,int seedy,
                      int color,int target,
                      double sigma,
                      Perlin& perlin,
                      std::mt19937& rng,
                      double warp_strength,
                      int fBm_octaves,
                      double persistence,
                      double lacunarity,
                      double freq_jitter,
                      double compactness_lambda){
    if(target<=0 || !inside(seedx,seedy,W,H) || img[seedy][seedx]!=-1) return;

    std::lognormal_distribution<double> logn(std::log(sigma),freq_jitter);
    double base_freq = 0.002 * logn(rng);

    double seedNoise = fBmWarp(perlin, seedx, seedy, base_freq,
                               fBm_octaves,persistence,lacunarity,warp_strength);

    std::priority_queue<Node,std::vector<Node>,NodeCmp> pq;
    pq.push({seedx,seedy,0.0});
    int placed=0;

    while(placed<target && !pq.empty()){
        auto [x,y,d]=pq.top(); pq.pop();
        if(!inside(x,y,W,H) || img[y][x]!=-1) continue;
        img[y][x]=color;
        placed++;
        for(int k=0;k<4;k++){
            int nx=x+DX4[k], ny=y+DY4[k];
            if(inside(nx,ny,W,H) && img[ny][nx]==-1){
                double nv=fBmWarp(perlin,nx,ny,base_freq,
                                  fBm_octaves,persistence,lacunarity,warp_strength);
                double dd=fabs(nv-seedNoise);
                if(compactness_lambda>0){
                    double dist = std::hypot(nx-seedx, ny-seedy);
                    dd += compactness_lambda * dist;
                }
                pq.push({nx,ny,dd});
            }
        }
    }
}

// ---------------------
// Main
// ---------------------
int main(int argc,char** argv){
    if(argc < 12){
        std::cerr<<"Usage: "<<argv[0]
                 <<" out.bmp W H alpha beta sigma p0 p1 p2 p3 seed min_patch"
                 <<" [warp_strength fBm_octaves persistence lacunarity freq_jitter compactness_lambda]\n";
        return 1;
    }
    std::string out=argv[1];
    int W=std::atoi(argv[2]);
    int H=std::atoi(argv[3]);
    double alpha=std::atof(argv[4]);
    double beta =std::atof(argv[5]); // unused placeholder
    double sigma=std::atof(argv[6]);
    double p0=std::atof(argv[7]);
    double p1=std::atof(argv[8]);
    double p2=std::atof(argv[9]);
    double p3=std::atof(argv[10]);
    unsigned int seed=std::atoi(argv[11]);
    int min_patch=std::atoi(argv[12]);

    // new parameters with defaults
    double warp_strength   = (argc>13)?std::atof(argv[13]):30.0;
    int fBm_octaves        = (argc>14)?std::atoi(argv[14]):4;
    double persistence     = (argc>15)?std::atof(argv[15]):0.5;
    double lacunarity      = (argc>16)?std::atof(argv[16]):2.0;
    double freq_jitter     = (argc>17)?std::atof(argv[17]):0.3;
    double compactness_lambda = (argc>18)?std::atof(argv[18]):0.0;

    std::cerr<<"Generating "<<W<<"x"<<H
             <<" alpha="<<alpha
             <<" sigma="<<sigma
             <<" min_patch="<<min_patch
             <<" warp="<<warp_strength
             <<" oct="<<fBm_octaves
             <<" pers="<<persistence
             <<" lac="<<lacunarity
             <<" fjit="<<freq_jitter
             <<" comp="<<compactness_lambda<<"\n";

    std::array<double,4> P={p0,p1,p2,p3};
    double sP=P[0]+P[1]+P[2]+P[3];
    if(sP<=0){P={0.25,0.25,0.25,0.25};sP=1.0;}
    for(int c=0;c<4;c++) P[c]/=sP;

    std::mt19937 rng(seed);
    Perlin perlin(seed);

    const int TOTAL=W*H;
    std::array<int,4> quota{};
    int acc=0;
    for(int c=0;c<4;c++){ quota[c]=(int)std::llround(P[c]*TOTAL); acc+=quota[c]; }
    while(acc<TOTAL){ quota[acc%4]++; acc++; }
    while(acc>TOTAL){ if(quota[acc%4]>0) quota[acc%4]--; acc--; }

    double tau=alpha_to_tau(alpha);
    std::uniform_int_distribution<int> RX(0,W-1), RY(0,H-1);

    struct PatchSpec{int color; int size;};
    std::vector<PatchSpec> patches;

    for(int c=0;c<4;c++){
        int remaining=quota[c];
        int hard_cap=std::max(min_patch,TOTAL/2);
        while(remaining>=min_patch){
            int s=sample_pareto_size(rng,min_patch,tau,hard_cap);
            if(s>remaining) s=remaining;
            patches.push_back({c,s});
            remaining-=s;
        }
        if(remaining>0){
            int s=std::max(min_patch,remaining);
            s=std::min(s,quota[c]);
            patches.push_back({c,s});
            remaining=0;
        }
    }
    std::sort(patches.begin(),patches.end(),
              [](const PatchSpec&a,const PatchSpec&b){return a.size>b.size;});

    std::vector<std::vector<int>> img(H,std::vector<int>(W,-1));
    for(auto &spec:patches){
        for(int attempt=0;attempt<200;++attempt){
            int sx=RX(rng), sy=RY(rng);
            if(img[sy][sx]!=-1) continue;
            grow_blob_perlin(img,W,H,sx,sy,spec.color,spec.size,
                             sigma,perlin,rng,
                             warp_strength,fBm_octaves,
                             persistence,lacunarity,
                             freq_jitter,compactness_lambda);
            break;
        }
    }

    std::discrete_distribution<int> pick({P[0],P[1],P[2],P[3]});
    for(int y=0;y<H;y++)
        for(int x=0;x<W;x++)
            if(img[y][x]==-1) img[y][x]=pick(rng);

    enforce_min_patch(img,W,H,min_patch);
    writeBMP(out,img,W,H);
    return 0;
}
