#!/usr/bin/env python3
"""verify_ln4_remarks.py — checks for the note on Entries 12.2.3 and 19.2.2 of Ramanujan's Lost Notebook, Part IV
(Andrews–Berndt, Springer 2013, pp. 270 and 377–380).   Requires numpy and mpmath.   Exit code 0 = all checks pass.

Part A (Entry 12.2.3, r = 0):  chi = the non-principal character mod 4,  A(x) = sum_{n<=x} chi(n) d(n).
Part B (Entry 19.2.2):  S(alpha) = (1/2) sum_d d^s cot(d alpha/2), the Abel value of sum sigma_s(n) sin(n alpha), s < -2.
"""
import math
import sys

import mpmath as mp
import numpy as np

RESULTS = []
rng = np.random.default_rng(20260927)  # same seed and draw order as the original computation (isbotlar/29)


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  —  {detail}" if detail else ""))


# ---------------- Part A ----------------
def S4(m):
    r = m % 4
    return ((r == 1) | (r == 2)).astype(np.int64)


def A(x):  # exact, O(sqrt x): chi(n)d(n) = sum_{ab=n} chi(a)chi(b), hyperbola method
    x = int(x)
    k = math.isqrt(x)
    a = np.arange(1, k + 1, dtype=np.int64)
    chi = np.where(a % 2 == 0, 0, np.where(a % 4 == 1, 1, -1))
    return int(2 * np.sum(chi * S4(x // a)) - int(S4(np.array([k]))[0]) ** 2)


mp.mp.dps = 30
chi4 = [0, 1, 0, -1]
L0 = mp.dirichlet(0, chi4)
check("A1 L(0,chi) = 1/2, so the residue of L(s,chi)^2 y^s/s at s=0 is 1/4", abs(L0 - mp.mpf(1) / 2) < mp.mpf(10) ** -25, mp.nstr(L0, 20))

Lam = lambda s: (mp.mpf(4) / mp.pi) ** s * mp.gamma((s + 1) / 2) ** 2 * mp.dirichlet(s, chi4) ** 2  # noqa: E731
pts = [mp.mpf("0.3"), mp.mpc("0.2", "3.1"), mp.mpc("-0.7", "1.4"), mp.mpf("1.9")]
worst = max(abs(Lam(s) / Lam(1 - mp.conj(s)).conjugate() - 1) if mp.im(s) else abs(Lam(s) / Lam(1 - s) - 1) for s in pts)
check("A2 Lambda(s) = (4/pi)^s Gamma((s+1)/2)^2 L(s,chi)^2 satisfies Lambda(s) = Lambda(1-s) (degree 2, omega = 1)",
      worst < mp.mpf(10) ** -20, f"max rel. diff {mp.nstr(worst, 3)}")

N = 10**6
d = np.zeros(N + 1, dtype=np.int64)
for i in range(1, N + 1):
    d[i::i] += 1
n = np.arange(N + 1)
chiN = np.where(n % 2 == 0, 0, np.where(n % 4 == 1, 1, -1))
Acum = np.cumsum(chiN * d)
xs = rng.integers(1, N, 400)
check("A3 hyperbola A(x) equals the sieve value (400 random x <= 10^6)", all(A(x) == Acum[x] for x in xs))

means = []
for X in (10**6, 10**8, 10**10, 10**12):
    xs = rng.integers(X // 2, X, 600)
    v = np.array([A(x) for x in xs], dtype=float)
    means.append(float(np.mean(v ** 2 / np.sqrt(xs))))
check("A4 mean of A(x)^2/sqrt(x) over 600 random x in [X/2,X], X = 10^6..10^12, equals the values in Table 1",
      [round(m, 4) for m in means] == [0.7322, 0.7815, 0.7571, 0.8151], ", ".join(f"{m:.4f}" for m in means))

T = np.cumsum(chiN[1:] * d[1:] * n[1:].astype(float) ** (-0.2))
osc = [float(T[a - 1:b].max() - T[a - 1:b].min()) for a, b in ((10**3, 10**4), (10**4, 10**5), (10**5, 10**6))]
check("A5 partial sums of sum chi(n)d(n)n^{-1/5}: oscillation over [10^k,10^{k+1}] grows (k = 3,4,5)",
      osc[0] < osc[1] < osc[2], ", ".join(f"{o:.3f}" for o in osc))


means_r = []
for rr in (-0.25, 0.25):
    Nr = 10**6
    sg = np.zeros(Nr + 1)
    for dv in range(1, Nr + 1):
        sg[dv::dv] += dv ** rr
    Ac = np.cumsum(chiN * sg)
    for X in (10**5, 10**6):
        xs = rng.integers(X // 2, X, 4000)
        means_r.append(float(np.mean(Ac[xs] ** 2 / xs ** (0.5 + rr))))
check("A6 (Remark 2, not a theorem) mean of A_r(x)^2/x^{1/2+r} is of order 1 for r = -1/4, 1/4 and X = 10^5, 10^6",
      all(0.5 < m < 1.5 for m in means_r), ", ".join(f"{m:.4f}" for m in means_r))

# ---------------- Part B ----------------
def S_sum(a, s, D):
    dd = np.arange(1, D + 1, dtype=np.float64)
    return 0.5 * float(np.sum(dd ** s / np.tan(np.mod(dd * a / 2, np.pi))))


def Lside(a, s, Sa):
    return a ** ((s + 1) / 2) * (float(mp.zeta(1 - s)) / a + 0.5 * float(mp.zeta(-s)) * math.tan(math.pi * s / 2) + Sa)


def min_dnorm(xi, D):  # min_{d<=D} d*||d xi||  (lower constant c in ||d xi|| >= c/d)
    dd = np.arange(1, D + 1, dtype=np.float64)
    fr = np.abs(dd * xi - np.round(dd * xi))
    return float(np.min(dd * fr))


s, xi = -2.5, math.sqrt(2)
a = 2 * math.pi * xi
Nn = 200000
sig = np.zeros(Nn + 1)
for dv in range(1, Nn + 1):
    sig[dv::dv] += dv ** s
r = 1 - 2e-4
nn = np.arange(1, Nn + 1)
abel = float(np.sum(sig[1:] * np.sin(nn * a) * r ** nn))
S_ref = S_sum(a, s, 4_000_000)
check("B1 Abel sum sum sigma_s(n) sin(n alpha) r^n (r = 1 - 2e-4) is close to S(alpha)  (s = -2.5, alpha/2pi = sqrt 2)",
      abs(abel - S_ref) < 1e-3, f"Abel {abel:.8f}  S {S_ref:.8f}")

D = 4_000_000
rows, ok = [], True
for xi in (math.sqrt(2), (1 + math.sqrt(5)) / 2, math.sqrt(3)):
    for s in (-2.5, -4.5):
        a = 2 * math.pi * xi
        b = 4 * math.pi ** 2 / a
        ca = cb = 0.25  # ||d xi|| > 1/((M+2)d), partial quotients of xi and 1/xi bounded by M <= 2 (Khinchin)
        # tail of (1/2) sum_{d>D} d^s |cot(pi d xi)| <= sum_{d>D} d^{s+1}/(4c) <= D^{s+2}/((-s-2) 4c)
        tail = lambda c: D ** (s + 2) / ((-s - 2) * 4 * c)  # noqa: E731
        diff = Lside(a, s, S_sum(a, s, D)) - Lside(b, s, S_sum(b, s, D))
        err = a ** ((s + 1) / 2) * tail(ca) + b ** ((s + 1) / 2) * tail(cb)
        rows.append(f"xi={xi:.4f} s={s}: diff={diff:+.4f} (tail<{err:.1e})")
        ok &= abs(diff) > 10 * err
check("B2 (19.2.2) fails for the Abel values: |L(alpha) - L(beta)| exceeds 10x the tail bound (3 quadratic irrationals x 2 values of s)",
      ok, "\n      " + "\n      ".join(rows))
check("B3 the constant c = 1/4 used in B2 is consistent: min_{d<=10^6} d*||d xi|| >= 1/4 for all six numbers",
      all(min_dnorm(x, 10**6) >= 0.25 for x in (math.sqrt(2), 1 / math.sqrt(2), (1 + math.sqrt(5)) / 2, 2 / (1 + math.sqrt(5)), math.sqrt(3), 1 / math.sqrt(3))))

s = -2.5
z = float(mp.zeta(1 - s))
rel = []
for p, q in ((1, 2), (1, 3), (2, 5), (3, 7)):
    a0 = 2 * math.pi * p / q
    e = 1e-6
    rel.append(abs(e * S_sum(a0 + e, s, 4_000_000) / (q ** (s - 1) * z) - 1))
check("B4 near alpha0 = 2 pi p/q: (alpha - alpha0) S(alpha) ~ q^{s-1} zeta(1-s)  (p/q = 1/2, 1/3, 2/5, 3/7; eps = 1e-6)",
      max(rel) < 1e-3, "rel. diffs " + ", ".join(f"{x:.1e}" for x in rel))

print(f"\n{sum(RESULTS)}/{len(RESULTS)} checks passed")
sys.exit(0 if all(RESULTS) else 1)
