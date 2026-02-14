#!/bin/bash
#SBATCH --account=2016345
#SBATCH --partition=biocomp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=03:00:00
#SBATCH --job-name=spa_mosaic
#SBATCH --mail-user=mfricke@unm.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=24GB
#SBATCH --output=debug_%j.out
#SBATCH --error=debug_%j.out

module load miniconda3
source activate substrate_generator

cd $SLURM_SUBMIT_DIR

OUT=out_spa_512
mkdir -p $OUT/thumbs

# --- Generate thumbnails ---
python3 - <<'EOF'
import os
from PIL import Image, ImageDraw, ImageFont

OUT="out_spa_512"
src_dir=f"{OUT}/bmps"
dst_dir=f"{OUT}/thumbs"
os.makedirs(dst_dir, exist_ok=True)
W,H=800,800
FONT_CANDIDATES=[
    "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/Library/Fonts/Menlo.ttc",
    "/System/Library/Fonts/SFNSMono.ttf"
]
font_path=next((p for p in FONT_CANDIDATES if os.path.exists(p)),None)
def load_font(px):
    from PIL import ImageFont
    if font_path:
        try: return ImageFont.truetype(font_path,px)
        except: pass
    return ImageFont.load_default()
for f in sorted(os.listdir(src_dir)):
    if not f.lower().endswith(".bmp"): continue
    src=os.path.join(src_dir,f)
    dst=os.path.join(dst_dir,f.replace(".bmp","_thumb.jpg"))
    with Image.open(src) as im:
        im.thumbnail((W,H))
        draw=ImageDraw.Draw(im,"RGBA")
        text=" ".join(p for p in f.split("_") if p.startswith(("a","s","m","w")))
        font=load_font(max(12,int(H*0.06)))
        bbox=font.getbbox(text); tw,th=bbox[2]-bbox[0],bbox[3]-bbox[1]
        draw.rectangle([0,H-th-12,tw+24,H],fill=(0,0,0,180))
        draw.text((12,H-th-8),text,font=font,fill=(255,255,255,255))
        im.save(dst,"JPEG",quality=90)
print("Thumbnails complete.")
EOF

# --- Assemble mosaic ---
python3 - <<'EOF'
import os
from PIL import Image, ImageDraw, ImageFont

OUT="out_spa_512"
thumb_dir=f"{OUT}/thumbs"
out_dir=OUT
files=sorted([os.path.join(thumb_dir,f) for f in os.listdir(thumb_dir) if f.endswith("_thumb.jpg")])
cols,rows=32,16
assert len(files)==cols*rows, f"Expected {cols*rows} thumbs, got {len(files)}"
imgs=[Image.open(f).convert('RGB') for f in files]
w,h=imgs[0].size
mosaic=Image.new('RGB',(cols*w+(cols-1)*4,rows*h+(rows-1)*4),'white')

FONT_CANDIDATES=[
    "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/Library/Fonts/Menlo.ttc",
    "/System/Library/Fonts/SFNSMono.ttf"
]
font_path=next((p for p in FONT_CANDIDATES if os.path.exists(p)),None)
def load_font(px):
    from PIL import ImageFont
    if font_path:
        try: return ImageFont.truetype(font_path,px)
        except: pass
    return ImageFont.load_default()
font=load_font(max(10,int(h*0.06)))
k=0
for r in range(rows):
    for c in range(cols):
        x,y=c*(w+4),r*(h+4)
        im=imgs[k]; mosaic.paste(im,(x,y))
        draw=ImageDraw.Draw(mosaic,"RGBA")
        text=f"[{k}]"
        bbox=font.getbbox(text); tw,th=bbox[2]-bbox[0],bbox[3]-bbox[1]
        tx,ty=x+(w-tw)//2,y+h-th-2
        draw.rectangle([tx-4,ty-2,tx+tw+4,ty+th+2],fill=(0,0,0,150))
        draw.text((tx,ty),text,font=font,fill=(255,255,255,255))
        k+=1
mosaic.save(os.path.join(out_dir,"spa_mosaic_all_512.jpg"),quality=95)
print("Mosaic written.")
EOF
