#!/usr/bin/env python3
"""Figure for the note on Entries 12.2.3 and 19.2.2.  Output: fig/ln4_remarks.png, fig/ln4_remarks.pdf"""
import math
import os

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42  # TrueType (arXiv flags Type 3 fonts)
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

os.makedirs("fig", exist_ok=True)

# left: partial sums T_s(x) = sum_{n<=x} chi(n) d(n) n^{-s}
N = 10**6
d = np.zeros(N + 1, dtype=np.int64)
for i in range(1, N + 1):
    d[i::i] += 1
n = np.arange(N + 1)
chi = np.where(n % 2 == 0, 0, np.where(n % 4 == 1, 1, -1))
T = {s: np.cumsum(chi[1:] * d[1:] * n[1:].astype(float) ** (-s)) for s in (0.2, 0.4)}
xs = np.unique(np.logspace(1, 6, 4000).astype(int)) - 1

# right: S(alpha) = (1/2) sum_d d^s cot(d alpha/2), s = -2.5, on xi in [0.30, 0.70]
s = -2.5
D = 20000
dd = np.arange(1, D + 1, dtype=np.float64)
xi = np.linspace(0.30, 0.70, 6001) + 1e-7 * math.sqrt(2)  # shifted off exact rationals
Sv = np.array([0.5 * np.sum(dd ** s / np.tan(np.mod(dd * math.pi * x, math.pi))) for x in xi])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
ax1.plot(xs + 1, T[0.2][xs], lw=0.7, color="#d62728", label=r"$s=0.2$ (diverges: $s<1/4$)")
ax1.plot(xs + 1, T[0.4][xs], lw=1.2, color="#1f77b4", label=r"$s=0.4$ (converges: $s>1/3$)")
ax1.set_xscale("log")
ax1.set_xlabel("x")
ax1.set_ylabel(r"$\sum_{n\leq x}\chi(n)d(n)n^{-s}$")
ax1.set_title("Entry 12.2.3, r = 0: partial sums", fontsize=10)
ax1.legend(fontsize=8.5, frameon=False, loc="upper left")

ax2.plot(xi, np.clip(Sv, -3, 3), lw=0.7, color="#333")
for p, q in ((1, 3), (2, 5), (3, 7), (1, 2), (4, 7), (3, 5), (2, 3)):
    if 0.30 <= p / q <= 0.70:
        ax2.axvline(p / q, color="#d62728", lw=0.6, ls=":")
        ax2.text(p / q, 2.75, f"{p}/{q}", color="#d62728", fontsize=8, ha="center")
ax2.set_ylim(-3, 3)
ax2.set_xlabel(r"$\alpha/2\pi$")
ax2.set_ylabel(r"$S(\alpha)$  (clipped to $[-3,3]$)")
ax2.set_title(r"Entry 19.2.2, $s=-2.5$: Abel value $S(\alpha)=\frac{1}{2}\sum_d d^s\cot(d\alpha/2)$", fontsize=10)

fig.tight_layout()
fig.savefig("fig/ln4_remarks.png", dpi=200)
fig.savefig("fig/ln4_remarks.pdf")
print("T_0.2 range on [1e5,1e6]:", round(float(T[0.2][10**5 - 1:].max() - T[0.2][10**5 - 1:].min()), 3),
      "| T_0.4 range on [1e5,1e6]:", round(float(T[0.4][10**5 - 1:].max() - T[0.4][10**5 - 1:].min()), 3))
