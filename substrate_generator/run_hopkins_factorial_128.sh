#!/usr/bin/env bash
set -euo pipefail

# Output folders
OUT=out_spa_512
mkdir -p "$OUT"/{bmps,thumbs,meta}

# Canvas + thumbs
W=5000
H=5000
THUMB_W=800
THUMB_H=800

# Discrete factors (balanced across 512 via stratified sampling)
ALPHAS=(0.6 0.8 1.0 1.2)                  # 4 values → each 128x
SIGMAS=(1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0)  # 8 values → each 64x
MIN_PATCHES=(100 500 2000 8000)           # 4 values → each 128x
PREVALS=("0.25 0.25 0.25 0.25")           # fixed balance

# “Extra” continuous params (explored via LHS ranges)
#   warp_strength:    0 .. 60 (px)
#   fBm_octaves:      3 .. 6 (integer)
#   persistence:      0.40 .. 0.70
#   lacunarity:       1.80 .. 2.40
#   freq_jitter:      0.00 .. 0.50
#   compactness:      0.000 .. 0.010
LHS_WARP_MIN=0
LHS_WARP_MAX=60
LHS_OCT_MIN=3
LHS_OCT_MAX=6
LHS_PERS_MIN=0.40
LHS_PERS_MAX=0.70
LHS_LAC_MIN=1.80
LHS_LAC_MAX=2.40
LHS_FJIT_MIN=0.00
LHS_FJIT_MAX=0.50
LHS_COMP_MIN=0.000
LHS_COMP_MAX=0.010

# Unused placeholder kept for compatibility
BETA=1.0

# Base seed; we’ll use SEED+idx for per-run determinism
SEED=12345

MANIFEST="$OUT/meta/manifest.csv"
echo "index,filename,width,height,alpha,beta,sigma,min_patch,warp,oct,pers,lac,fjit,comp,p0,p1,p2,p3,seed" > "$MANIFEST"

# --------------------------
# Step 0: generate a Latin Hypercube (512 rows) into a CSV the shell will read
# --------------------------
python3 - <<'EOF' > out_spa_512/meta/params_512.csv
import csv, random, math

N = 512
random.seed(123456)

# Discrete axes, balanced across N
ALPHAS = [0.6, 0.8, 1.0, 1.2]          # 4 → each 128x
SIGMAS = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0]  # 8 → each 64x
MINPATCH = [100,500,2000,8000]         # 4 → each 128x
P = (0.25,0.25,0.25,0.25)

# Build stratified lists
alphas = (ALPHAS * (N // len(ALPHAS)))
sigmas = (SIGMAS * (N // len(SIGMAS)))
minps  = (MINPATCH * (N // len(MINPATCH)))

random.shuffle(alphas)
random.shuffle(sigmas)
random.shuffle(minps)

# LHS for continuous params: generate N strata in [0,1), permute, then map to ranges
def lhs_strata(N):
    bins = [(i + random.random())/N for i in range(N)]
    random.shuffle(bins)
    return bins

u_warp = lhs_strata(N)
u_oct  = lhs_strata(N)
u_pers = lhs_strata(N)
u_lac  = lhs_strata(N)
u_fjit = lhs_strata(N)
u_comp = lhs_strata(N)

# Ranges
WARP_MIN, WARP_MAX = 0.0, 60.0
OCT_MIN, OCT_MAX   = 3, 6
PERS_MIN, PERS_MAX = 0.40, 0.70
LAC_MIN,  LAC_MAX  = 1.80, 2.40
FJIT_MIN, FJIT_MAX = 0.00, 0.50
COMP_MIN, COMP_MAX = 0.000, 0.010

rows = []
for i in range(N):
    warp = WARP_MIN + u_warp[i]*(WARP_MAX - WARP_MIN)
    octv = int(round(OCT_MIN + u_oct[i]*(OCT_MAX - OCT_MIN)))
    pers = PERS_MIN + u_pers[i]*(PERS_MAX - PERS_MIN)
    lac  = LAC_MIN  + u_lac[i]*(LAC_MAX  - LAC_MIN)
    fjit = FJIT_MIN + u_fjit[i]*(FJIT_MAX - FJIT_MIN)
    comp = COMP_MIN + u_comp[i]*(COMP_MAX - COMP_MIN)

    rows.append({
        "alpha": alphas[i],
        "sigma": sigmas[i],
        "min_patch": minps[i],
        "warp": warp,
        "oct": octv,
        "pers": pers,
        "lac": lac,
        "fjit": fjit,
        "comp": comp
    })

# Stable ordering: by min_patch, then alpha, then sigma (so mosaics read as a map)
rows.sort(key=lambda r: (r["min_patch"], r["alpha"], r["sigma"], r["warp"]))

with open("out_spa_512/meta/params_512.csv","w",newline="") as f:
    w = csv.writer(f)
    w.writerow(["alpha","sigma","min_patch","warp","oct","pers","lac","fjit","comp"])
    for r in rows:
        w.writerow([r["alpha"],r["sigma"],r["min_patch"],r["warp"],r["oct"],r["pers"],r["lac"],r["fjit"],r["comp"]])
EOF

# --------------------------
# Step 1: generate BMPs
# --------------------------
idx=0
while IFS=, read -r alpha sigma mp warp oct pers lac fjit comp; do
  # skip header
  if [[ "$alpha" == "alpha" ]]; then continue; fi

  read -r p0 p1 p2 p3 <<< "${PREVALS[0]}"

  # integerize/format where needed
  alpha_fmt=$(printf "%.2f" "$alpha")
  sigma_fmt=$(printf "%.1f" "$sigma")
  mp_fmt=$(printf "%05d" "$mp")
  warp_int=$(printf "%.0f" "$warp")
  oct_int=$oct

  fname=$(printf "spa_%03d_a%s_s%s_m%s_w%02d.bmp" \
                  "$idx" "$alpha_fmt" "$sigma_fmt" "$mp_fmt" "$warp_int")
  outpath="$OUT/bmps/$fname"

  if [[ ! -f "$outpath" ]]; then
    echo "[$idx] Generating $outpath"
    run_seed=$((SEED + idx))
    ./reef_spa_substrate "$outpath" $W $H \
        "$alpha" $BETA "$sigma" \
        $p0 $p1 $p2 $p3 $run_seed $mp \
        "$warp" "$oct" "$pers" "$lac" "$fjit" "$comp"
  fi

  echo "$idx,$fname,$W,$H,$alpha,$BETA,$sigma,$mp,$warp,$oct,$pers,$lac,$fjit,$comp,$p0,$p1,$p2,$p3,$((SEED+idx))" >> "$MANIFEST"
  idx=$((idx+1))
done < "$OUT/meta/params_512.csv"

echo "Generated $idx BMPs in $OUT/bmps"

# --------------------------
# Step 2: generate annotated thumbnails (scalable font, explicit path)
# --------------------------
python3 - <<'EOF'
import os
from PIL import Image, ImageDraw, ImageFont

OUT="out_spa_512"
src_dir=f"{OUT}/bmps"
dst_dir=f"{OUT}/thumbs"
os.makedirs(dst_dir, exist_ok=True)

W, H = 800, 800  # thumbnails
# Try explicit font paths first, then fallbacks
FONT_CANDIDATES = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",      # Linux common
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/Library/Fonts/Menlo.ttc",                                 # macOS
    "/System/Library/Fonts/SFNSMono.ttf",
    "DejaVuSansMono.ttf"
]
font_path = None
for cand in FONT_CANDIDATES:
    if os.path.exists(cand):
        font_path = cand
        break

def load_font(px):
    if font_path:
        try:
            return ImageFont.truetype(font_path, px)
        except Exception:
            pass
    return ImageFont.load_default()

for f in sorted(os.listdir(src_dir)):
    if not f.lower().endswith(".bmp"):
        continue
    src = os.path.join(src_dir, f)
    dst = os.path.join(dst_dir, f.replace(".bmp","_thumb.jpg"))
    with Image.open(src) as im:
        im.thumbnail((W,H))
        draw = ImageDraw.Draw(im, "RGBA")
        base = os.path.basename(f).replace(".bmp","")
        # Label with a.. s.. m.. (and warp as small hint)
        parts = [p for p in base.split("_") if p.startswith(("a","s","m","w"))]
        text = " ".join(parts)
        fs = max(12, int(H * 0.06))  # ~6% of height
        font = load_font(fs)
        bbox = font.getbbox(text)
        tw, th = bbox[2]-bbox[0], bbox[3]-bbox[1]
        draw.rectangle([0, H-th-12, tw+24, H], fill=(0,0,0,180))
        draw.text((12, H-th-8), text, font=font, fill=(255,255,255,255))
        im.save(dst, "JPEG", quality=90)
print("Annotated thumbnails written to", dst)
EOF

# --------------------------
# Step 3: assemble single mosaic (ALL 512 runs) in OUT/ (above bmps/)
# --------------------------
python3 - <<'EOF'
import os
from PIL import Image, ImageDraw, ImageFont

OUT="out_spa_512"
thumb_dir=f"{OUT}/thumbs"
out_dir=OUT  # one level above bmps

files = sorted([os.path.join(thumb_dir, f)
                for f in os.listdir(thumb_dir) if f.endswith("_thumb.jpg")])

# 32 x 16 = 512 panels
cols, rows = 32, 16
assert len(files) == cols*rows, f"Expected {cols*rows} thumbs, got {len(files)}"

imgs = [Image.open(f).convert('RGB') for f in files]
w,h = imgs[0].size
mosaic = Image.new('RGB', (cols*w+(cols-1)*4, rows*h+(rows-1)*4), 'white')

# Scalable font for labels (same robust path logic)
FONT_CANDIDATES = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/Library/Fonts/Menlo.ttc",
    "/System/Library/Fonts/SFNSMono.ttf",
    "DejaVuSansMono.ttf"
]
font_path = next((p for p in FONT_CANDIDATES if os.path.exists(p)), None)
def load_font(px):
    from PIL import ImageFont
    if font_path:
        try:
            return ImageFont.truetype(font_path, px)
        except Exception:
            pass
    return ImageFont.load_default()

fs = max(10, int(h * 0.06))
font = load_font(fs)

k=0
for r in range(rows):
    for c in range(cols):
        x,y = c*(w+4), r*(h+4)
        im = imgs[k]
        mosaic.paste(im, (x,y))

        # pull minimal label from filename already burned into thumb if you prefer;
        # but we'll add a light header showing index to help debugging
        draw = ImageDraw.Draw(mosaic, "RGBA")
        label = os.path.basename(files[k]).replace("_thumb.jpg","")
        # shorten label to key fields:
        parts = [p for p in label.split("_") if p.startswith(("a","s","m"))]
        text = " ".join(parts)
        bbox = font.getbbox(text)
        tw, th = bbox[2]-bbox[0], bbox[3]-bbox[1]
        tx, ty = x + (w - tw)//2, y + h - th - 2
        draw.rectangle([tx-4, ty-2, tx+tw+4, ty+th+2], fill=(0,0,0,150))
        draw.text((tx, ty), text, font=font, fill=(255,255,255,255))

        k+=1

out_path = os.path.join(out_dir, "spa_mosaic_all_512.jpg")
mosaic.save(out_path, quality=95)
print("Mosaic written to", out_path)
EOF

echo "Done."