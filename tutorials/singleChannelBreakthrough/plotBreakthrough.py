#!/usr/bin/env python3
"""Plot the outlet breakthrough curve and the mean loading and temperature of the wall.
Run after the simulation: python3 plotBreakthrough.py  (needs numpy and matplotlib)."""
import glob
import numpy as np
import matplotlib.pyplot as plt

YIN = 6.285e-4   # feed CO2 mass fraction, as in 0/fluid/CO2


def load(pattern):
    return np.loadtxt(sorted(glob.glob(pattern))[-1], comments="#")


b = load("postProcessing/fluid/outletCO2/*/surfaceFieldValue.dat")
s = load("postProcessing/solid/solidAverage/*/volFieldValue.dat")
fig, ax = plt.subplots(1, 3, figsize=(12, 3.5))
ax[0].plot(b[:, 0], b[:, 1]/YIN); ax[0].set_ylabel("outlet c/c0")
ax[1].plot(s[:, 0], s[:, 1]); ax[1].set_ylabel("mean loading q [mol per kg wall]")
ax[2].plot(s[:, 0], s[:, 2]); ax[2].set_ylabel("mean wall temperature [K]")
for a in ax:
    a.set_xlabel("time [s]")
fig.tight_layout(); fig.savefig("breakthrough.png", dpi=150)
print("written breakthrough.png")
