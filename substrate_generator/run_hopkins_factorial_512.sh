#!/usr/bin/env bash
set -euo pipefail

# ---------- config ----------
EXEC=./reef_multifield           # <- your compiled generator
OUT=out_hopkins_1024
mkdir -p "$OUT"/{bmps,thumbs,meta}

# Canvas: 5000 px ≈ 500 m → 0.1 m/px
W=5000
H=5000

# Thumbnails & mosaic
THUMB_W=420
THUMB_H=420
GAP=4
COLS=32
ROWS=32  # 32x32 = 1024

SEED=12345

# ---------- factorial: 8 x 8 x 4 x 4 = 1024 ----------
# MIN_PATCH in pixels (≈ m^2 by *0.01). 20 px ≈ 0.2 m² ... 10k px ≈ 100 m²
MIN_PATCH=(20 100 200 500 1000 2000 5000 10000)

# HOPKINS as contrast knob (used as pow(val, 1/hopkins))
# Lower ~ more contrast/clumpy; higher ~ flatter/uniform.
HOPKINS=(0.40 0.60 0.80 1.00 1.25 1.50 2.00 2.50)

# CORR_LEN: number of “noise cycles” across 500 m.
# Feature size ≈ 500 m / corr_len → {~2m, 5m, 10m, 25m}
CORR_LEN=(250 100 50 20)

# PREVALS: 4 skews so each ecological class dominates once (no fixed order bias).
# Order (fixed colors): [Yellow=Complex+Food, Blue=Simple+NoFood, Red=Complex+NoFood, Green=Simple+Food]
PREVALS=(
  "0.55 0.15 0.15 0.15"  # Yellow-dominant
  "0.15 0.55 0.15 0.15"  # Blue-dominant
  "0.15 0.15 0.55 0.15"  # Red-dominant
  "0.15 0.15 0.15 0.55"  # Green-dominant
)

# ---------- manifest ----------
MANIFEST="$OUT/meta/manifest.csv"
echo "index,filename,width,height,min_patch,hopkins,corr_len,p0,p1,p2,p3,seed" > "$MANIFEST"

# ---------- generate ----------
idx=0
for mp in "${MIN_PATCH[@]}"; do
  for h in "${HOPKINS[@]}"; do
    for cl in "${CORR_LEN[@]}"; do
      for p in "${PREVALS[@]}"; do
        read -r p0 p1 p2 p3 <<< "$p"
        fname=$(printf "reef_%04d_mp%05d_h%.2f_cl%03d_p%s-%s-%s-%s_seed%u.bmp" \
                        "$idx" "$mp" "$h" "$cl" "$p0" "$p1" "$p2" "$p3" "$SEED")
        outpath="$OUT/bmps/$fname"

        if [[ ! -f "$outpath" ]]; then
          echo "[$idx/1024] $outpath"
          "$EXEC" "$outpath" $W $H "$mp" "$h" "$cl" $p0 $p1 $p2 $p3 $SEED
        fi

        echo "$idx,$fname,$W,$H,$mp,$h,$cl,$p0,$p1,$p2,$p3,$SEED" >> "$MANIFEST"
        idx=$((idx+1))
      done
    done
  done
done
echo "Generated $idx BMPs in $OUT/bmps"

# ---------- thumbnails with labels ----------
export OUT THUMB_W THUMB_H

python3 - <<'PY'
import os
from PIL import Image, ImageDraw, ImageFont

OUT = os.environ["OUT"]
THUMB_W = int(os.environ["THUMB_W"])
THUMB_H = int(os.environ["THUMB_H"])
src_dir = os.path.join(OUT, "bmps")
dst_dir = os.path.join(OUT, "thumbs")
os.makedirs(dst_dir, exist_ok=True)

# try a readable font
font = None
for fp in [
    "/System/Library/Fonts/Supplemental/Courier New Bold.ttf",
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/Library/Fonts/Courier New.ttf",
    "DejaVuSansMono.ttf",
]:
    try:
        font = ImageFont.truetype(fp, 24)
        break
    except Exception:
        pass
if font is None:
    font = ImageFont.load_default()

files = sorted([f for f in os.listdir(src_dir) if f.lower().endswith(".bmp")])
for f in files:
    src = os.path.join(src_dir, f)
    dst = os.path.join(dst_dir, f.replace(".bmp","_thumb.jpg"))
    with Image.open(src) as im:
        im.thumbnail((THUMB_W, THUMB_H))
        draw = ImageDraw.Draw(im, "RGBA")
        text = f[:-4]
        bbox = draw.textbbox((0,0), text, font=font)
        tw, th = bbox[2]-bbox[0], bbox[3]-bbox[1]
        pad = 8
        draw.rectangle([0, THUMB_H-th-pad, min(tw+2*pad, THUMB_W), THUMB_H],
                       fill=(0,0,0,175))
        draw.text((pad, THUMB_H-th-pad//2), text, font=font,
                  fill=(255,255,255,255))
        im.save(dst, "JPEG", quality=90)
print("Annotated thumbnails ->", dst_dir)
PY

# ---------- mosaic (32 x 32 = 1024) ----------
export GAP COLS ROWS OUT

python3 - <<'PY'
import os
from PIL import Image

OUT = os.environ["OUT"]
GAP = int(os.environ.get("GAP","4"))
COLS = int(os.environ.get("COLS","32"))
ROWS = int(os.environ.get("ROWS","32"))
thumb_dir = os.path.join(OUT, "thumbs")
out_path  = os.path.join(OUT, "mosaic_32x32.jpg")

files = sorted([os.path.join(thumb_dir, f)
                for f in os.listdir(thumb_dir)
                if f.endswith("_thumb.jpg")])
n = len(files)
assert n == COLS*ROWS, f"Expected {COLS*ROWS} thumbs, found {n}"

imgs = [Image.open(f).convert('RGB') for f in files]
w,h = imgs[0].size
mosaic = Image.new('RGB',
                   (COLS*w+(COLS-1)*GAP, ROWS*h+(ROWS-1)*GAP),
                   'white')

k=0
for r in range(ROWS):
    for c in range(COLS):
        mosaic.paste(imgs[k], (c*(w+GAP), r*(h+GAP)))
        k+=1

mosaic.save(out_path, quality=95)
print("Mosaic ->", out_path)
PY