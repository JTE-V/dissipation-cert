---
title: 'dissipation-cert: a dependency-free reference implementation with analytic branch bounds for 1/(sin x - x) and its regularisation'
tags:
  - Python
  - numerical analysis
  - validated numerics
  - error bounds
  - special functions
authors:
  - name: JTE
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 8 November 2026
bibliography: paper.bib
---

# Summary

The function F(x)=1/(sinx-x) has a third-order pole at the origin, since `sin x - x ~ -x^3/6` and therefore
`F ~ -6/x^3`. Evaluating it directly near `x = 0` is dominated by cancellation.
Subtracting the two leading divergences — the exact rational terms `6/x^3` and
`(3/10)/x` — leaves R(x) = F(x) + 6/x^3 + (3/10)/x
= -11x/1400 - 17x^3/126000 - ... 
which is bounded and can be evaluated with a stated bound. `dissipation-cert`
implements this reduction, together with a companion diagnostic that reports,
across a grid, which branch's error term governs the result.

The package uses only the Python standard library (`fractions`). All arithmetic
on the series branch is exact rational arithmetic, so the branch bound can be
checked against the source by hand: the bound is the first dropped term,
`|H_K| x^(2K-3)`, times a tail factor `1/(1-x^2)` which is at most `4/3` on the
domain where the branch is used. Because that bound contains no working
precision, it is independent of the precision at which the rest of the
computation runs.

# Statement of need

Verification is part of research, not a preliminary to it: a computed value that
carries no account of its own error is a claim, not a result. When a quantity is
defined by a function with a pole, the account of the error is the part that
takes work, and it is the part that a reader has to be able to re-derive
independently.

For a large class of functions, general-purpose arbitrary-precision libraries
answer this need well. They attach a radius to every intermediate quantity and
propagate it, which is the right tool for an arbitrary expression. But the
radius they report is a consequence of the working precision the user selected.
If a series truncation is what governs the error, that fact is not returned as a
separately derivable quantity that can be written down, checked, and compared
against a bound obtained elsewhere.

`dissipation-cert` addresses that gap for one specific reduction. It separates
the two regimes of `R(x)` explicitly, and on the regime where the series
converges it returns a bound that is analytic, is stated as a single named
quantity, and does not depend on the working precision. The research purpose is
therefore twofold:

- to make one reduction available in a form whose bound can be checked without
  re-running anything, and can be compared against an independently written
  implementation — the CI workflow performs that comparison on a machine that is
  not the author's, for three Python versions;
- to serve as a worked example of a branch-and-bound evaluation in which the
  switch between branches, and its cost, is reported rather than left implicit.

The package is not a replacement for a general arbitrary-precision library, and
no head-to-head accuracy comparison against one has been performed.
# Comparison with existing software

`Arb` (now part of `FLINT`), `MPFI`, `INTLAB` and the interval-arithmetic
subsets of `JuliaIntervals` all provide rigorous enclosures of real functions,
including special functions, and all of them generalise well beyond the single
quantity considered here [@johansson2017arb; @flint; @mpfi; @rump1999intlab;
@ieee1788]. They track the error alongside the value, which makes them suitable
for arbitrary expressions; the cost is that the reported bound follows from the
working precision, and that the bound for any particular series truncation is
not visible as a separate, derivable quantity.

This package makes the opposite trade. It handles one reduction, hard-codes the
two branches, and states the series bound analytically and independently of
precision. In exchange it cannot handle anything but the function named in the
title.

No head-to-head accuracy or performance comparison against any of these
libraries has been carried out. Such a comparison is a natural next step and is
explicitly not claimed here.

# Notes on the bounds

Two properties are worth stating because they are easy to misread.

The bound on the series branch is analytic and precision-independent. The
quantity reported on the direct branch is not a bound: it is the observed
difference between two runs at different working precisions, and it is
therefore a magnitude estimate. This is stated in the documentation and is
visible in the printed output, but it is the weakest claim in the package and
should be read as such.

There is a discontinuity of tens of orders of magnitude between the two
branches at `x = 1/2`. It belongs to the switch between branches, not to
`R(x)`, which is analytic across that point. The diagnostic exists mainly to
make that step, and its size, explicit.
