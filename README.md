# Two open points in Part IV of Ramanujan's Lost Notebook

A short note on Entries 12.2.3 and 19.2.2 of Andrews–Berndt, *Ramanujan's Lost Notebook, Part IV* (Springer 2013).

- **Entry 12.2.3 (r = 0).** The authors expect the identity L(s,χ)² = Σ χ(n)d(n)n^(−s) (χ mod 4) to hold for Re s > 0.
  It does not: the abscissa of convergence of the series lies in [1/4, 1/3] (Theorem 1). The lower bound is a consequence of
  classical Ω-results; the note derives it from a mean-square theorem of Cao–Tanigawa–Zhai whose hypotheses are checked.
- **Entry 19.2.2.** Read by Abel summation (s < −2, α/2π a quadratic irrational) the divergent series has the value
  ½Σ d^s cot(dα/2), but the identity still fails (Corollary 5); numerically this value has poles at rational α/2π.

No new methods are claimed.

- **Note:** [`paper/Kenjaev_LN4_Entries_12_2_3_and_19_2_2.pdf`](paper/Kenjaev_LN4_Entries_12_2_3_and_19_2_2.pdf) (5 pages)
- **Calculator:** https://selloriybrendi.github.io/ramanujan-lost-notebook-iv-remarks/ (offline, built-in self-test)
- **Verification:** `python3 verify/verify_ln4_remarks.py` → `10/10 checks passed` (NumPy, mpmath)
- **Figure:** `python3 verify/make_figure.py` (needs matplotlib)

## Use of generative AI
Claude Opus 5.5 (Anthropic, `claude-opus-5-5`) via Claude Code was used to locate the open questions, to find and check the
mean-square theorem, for the computations (NumPy, mpmath) and to help prepare the text, under the author's direction; every number
was verified independently. The author takes full responsibility.

Author: Otakhon U. Kenjaev, ORCID [0009-0009-3566-9285](https://orcid.org/0009-0009-3566-9285)
