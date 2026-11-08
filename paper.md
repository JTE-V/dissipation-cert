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

Evaluating a quotient of this kind near a pole is a routine task, and
arbitrary-precision libraries perform it well. What they generally do not do is
separate two questions that this package separates. A ball-arithmetic library
attaches an error bound to every intermediate quantity and propagates it. That
is more general than is needed for a single fixed reduction, and it couples the
reported bound to the working precision: the user selects a precision, and the
radius reported at the end follows from that choice. For a series truncation,
the resulting bound is not available as a separately derivable quantity that
can be written down and checked on its own.

`dissipation-cert` treats the two regimes of the same function separately, and
for the regime where the series converges it states a bound that holds
regardless of precision. The intended use is narrow:

- as an independent second implementation of one particular reduction, so that
  a value obtained elsewhere can be checked against something written from
  scratch — with the check running on a machine that is not the author's (the
  CI workflow does this on GitHub's runners, for three Python versions);
- as a worked example of the branch-and-bound pattern, in which the branch
  switch, and its consequence — a step of tens of orders of magnitude — is
  reported rather than left implicit.

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
