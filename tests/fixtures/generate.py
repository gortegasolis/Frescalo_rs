#!/usr/bin/env python3
"""Generate the synthetic test fixture used by tests/golden.rs.

The data are artificial: 100 sites on a 10 x 10 grid, 40 analysed species
recorded in three periods with uneven recording effort, and 60 well-recorded
training-set species. Species occupancy follows smooth gradients across the
grid, and a few species increase or decrease between periods, so every stage
of the pipeline (neighbourhoods, effort standardisation, benchmarks, time
factors) is exercised.

The generated files are committed; re-running this script (seed fixed)
reproduces them exactly. Run from this directory:  python3 generate.py
"""
import math
import random

rng = random.Random(20260928)

NX, NY = 10, 10
sites = []
for y in range(NY):
    for x in range(NX):
        sites.append((f"S{y * NX + x + 1:03d}", x, y))

# Projected coordinates (metres, 10 km grid)
with open("Samp_locations.txt", "w", newline="\n") as f:
    for s, x, y in sites:
        f.write(f"{s} {300000 + 10000 * x} {600000 + 10000 * y}\n")

# Geographic coordinates (longitude latitude): spans 35-66.5 degrees N and
# crosses the antimeridian, for the geodesic-mode fixture.
with open("Samp_lonlat.txt", "w", newline="\n") as f:
    for s, x, y in sites:
        lon = 175.0 + 1.5 * x
        if lon >= 180.0:
            lon -= 360.0
        lat = 35.0 + 3.5 * y
        f.write(f"{s} {lon:.4f} {lat:.4f}\n")


def logistic(v):
    return 1.0 / (1.0 + math.exp(-v))


def gradient_species(n):
    out = []
    for _ in range(n):
        out.append((rng.uniform(-2.5, 1.5), rng.uniform(-3, 3), rng.uniform(-3, 3)))
    return out


# Training set: well recorded, all periods pooled
training = gradient_species(60)
with open("Training.txt", "w", newline="\n") as f:
    for s, x, y in sites:
        for j, (a, b, c) in enumerate(training):
            p = logistic(a + b * x / (NX - 1) + c * y / (NY - 1))
            if rng.random() < p * 0.9:
                f.write(f"{s} Tr{j + 1:02d}\n")

# Analysed species, three periods
periods = ["P1990", "P2000", "P2010"]
species = gradient_species(40)
trend = [0.0] * 40
for j in range(0, 40, 7):
    trend[j] = 1.2      # increasing
for j in range(3, 40, 9):
    trend[j] = -1.2     # decreasing
detect = [rng.uniform(0.4, 1.0) for _ in species]
effort = {}
for s, _, _ in sites:
    for t in periods:
        e = rng.uniform(0.2, 1.0)
        if rng.random() < 0.1:
            e = 0.05            # barely visited
        effort[(s, t)] = e

with open("Test.txt", "w", newline="\n") as f:
    for s, x, y in sites:
        for ti, t in enumerate(periods):
            for j, (a, b, c) in enumerate(species):
                occ = logistic(a + b * x / (NX - 1) + c * y / (NY - 1) + trend[j] * (ti - 1))
                if rng.random() < occ * effort[(s, t)] * detect[j]:
                    f.write(f"{s} Sp{j + 1:02d} {t}\n")

with open("NotBench.txt", "w", newline="\n") as f:
    f.write("Sp05\nSp12\n")
