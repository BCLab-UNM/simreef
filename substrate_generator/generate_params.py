#!/usr/bin/env python3
"""
Generate Latin Hypercube–stratified parameter combinations
for the reef substrate simulation.

Output:
  out_spa_512/meta/params_512.csv
"""

import csv
import random
import os

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
N = 512                     # Number of total samples
SEED = 123456               # RNG seed for reproducibility

ALPHAS = [0.6, 0.8, 1.0, 1.2]
SIGMAS = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
MINPATCH = [100, 500, 2000, 8000]

# Continuous parameter ranges (for Latin Hypercube Sampling)
RANGES = {
    "warp": (0.0, 60.0),     # pixel warp strength
    "oct": (3, 6),           # fBm octaves (integer)
    "pers": (0.40, 0.70),    # fBm persistence
    "lac": (1.80, 2.40),     # fBm lacunarity
    "fjit": (0.00, 0.50),    # frequency jitter
    "comp": (0.000, 0.010)   # compactness lambda
}

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def lhs_strata(N: int):
    """Generate a Latin Hypercube vector in [0,1)."""
    bins = [(i + random.random()) / N for i in range(N)]
    random.shuffle(bins)
    return bins

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    random.seed(SEED)
    os.makedirs("out_spa_512/meta", exist_ok=True)

    # LHS uniform samples
    u_warp = lhs_strata(N)
    u_oct = lhs_strata(N)
    u_pers = lhs_strata(N)
    u_lac = lhs_strata(N)
    u_fjit = lhs_strata(N)
    u_comp = lhs_strata(N)

    # Expand and shuffle discrete factors so we have 512 combinations total
    alphas = (ALPHAS * (N // len(ALPHAS)))
    sigmas = (SIGMAS * (N // len(SIGMAS)))
    minps  = (MINPATCH * (N // len(MINPATCH)))

    random.shuffle(alphas)
    random.shuffle(sigmas)
    random.shuffle(minps)

    rows = []
    for i in range(N):
        warp = RANGES["warp"][0] + u_warp[i] * (RANGES["warp"][1] - RANGES["warp"][0])
        octv = int(round(RANGES["oct"][0] + u_oct[i] * (RANGES["oct"][1] - RANGES["oct"][0])))
        pers = RANGES["pers"][0] + u_pers[i] * (RANGES["pers"][1] - RANGES["pers"][0])
        lac  = RANGES["lac"][0]  + u_lac[i]  * (RANGES["lac"][1]  - RANGES["lac"][0])
        fjit = RANGES["fjit"][0] + u_fjit[i] * (RANGES["fjit"][1] - RANGES["fjit"][0])
        comp = RANGES["comp"][0] + u_comp[i] * (RANGES["comp"][1] - RANGES["comp"][0])

        rows.append([alphas[i], sigmas[i], minps[i], warp, octv, pers, lac, fjit, comp])

    # Sort for visually meaningful mosaic ordering (by min_patch, then alpha, etc.)
    rows.sort(key=lambda r: (r[2], r[0], r[1], r[3]))

    out_path = "out_spa_512/meta/params_512.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["alpha", "sigma", "min_patch", "warp", "oct", "pers", "lac", "fjit", "comp"])
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_path}")

# ---------------------------------------------------------------------
if __name__ == "__main__":
    main()
