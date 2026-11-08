# dissipation-cert

[![CI](https://github.com/JTE-V/dissipation-cert/actions/workflows/ci.yml/badge.svg)](https://github.com/JTE-V/dissipation-cert/actions/workflows/ci.yml)

A **standard-library-only** (Python 3, `fractions`) reference implementation that
evaluates

```
F(x)   = 1 / (sin x - x)
R(x)   = F(x) + 6/x^3 + (3/10)/x
Phi(x) = R(x) - Si(x)
```

together with **piecewise error bounds**, plus a diagnostic that reports which
branch's error term dominates at a given `x`.

## Statement of need

`F(x)` has a third-order pole at the origin: `sin x - x ~ -x^3/6`, hence
`F ~ -6/x^3`. Naive evaluation is therefore unusable near `x = 0`. Subtracting
the two leading divergences — the exact rational terms `6/x^3` and `(3/10)/x` —
leaves

```
R(x) = -11x/1400 - 17x^3/126000 - ...
```

which is bounded (`|R(x)| < 0.1` over the sampled range; `R(x) < 0` is proved
for `x >= 4`) and can be evaluated with a stated bound.

The implementation is deliberately dependency-free so that the series bound can
be audited against the source by hand. It is intended as a second, independent
implementation for auditing other evaluations of the same quantity — not as a
replacement for arbitrary-precision libraries.

## Requirements

- Python 3.11 or newer (CI-tested on 3.11, 3.12, 3.13)
- No third-party packages. Standard library `fractions` only.

## Usage

```
python3 dissipation_cert.py exact            # coefficients of F(x)
python3 dissipation_cert.py selfcheck        # internal consistency checks
python3 dissipation_cert.py bound            # piecewise bound report
python3 dissipation_cert.py probe 0.5        # value + bound at one point
python3 dissipation_cert.py scan 200000 30   # long sweep (minutes)

python3 stair_diag.py                        # error-vs-x diagnostic;
                                             # writes stair_diag.csv, stair_diag.svg
```

`stair_diag.py` imports `dissipation_cert`, so the two files must sit in the
same directory.

## Branch structure

| Range | Branch | Evaluation | Error term | Rigorous? |
|---|---|---|---|---|
| `x <= 1/2` | series | `R_series(K)`, retains `H_1 .. H_(K-1)` | first dropped term `\|H_K\| * x^(2K-3)` times a tail factor `1/(1-x^2)` | yes, analytic |
| `x > 1/2` | direct | evaluates `F`, `Si` directly | observed difference between two `prec` levels | no, magnitude estimate |

`NK = 24`, `SAFE = 1/2`. The series bound is **independent of `prec`**; the
direct branch's estimate is not.

## Expected output

Rows 1-5 and 7-8 are produced by the commands in column 2. Row 6 is an optional
external cross-check that requires `mpmath` and is **not** a dependency of this
package. The CI job reruns all in-repository commands on GitHub's runners and
uploads `stair_diag.csv` and `stair_diag.svg` as artifacts.

| # | Check | Command | Expected |
|---|---|---|---|
| 1 | power-series coefficients | `exact` | `-6`, `-3/10`, `-11/1400`, `-17/126000` |
| 2 | series bound at `x = 1/2` (tail factor included) | `selfcheck` | `4.3e-56` |
| 3 | first dropped term only at `x = 1/2` | `selfcheck` | `3.2e-56`  (= row 2 divided by `4/3`) |
| 4 | direct-branch collision vs `sin` first dropped term | `selfcheck` | ratio `0.987 ~ 0.998` |
| 5 | coefficient collision between the two implementations | `stair_diag` | all `27` coefficients equal |
| 6 | `actual / bound` at `x = 1/2` | optional, external `mpmath` 300 dps | `0.750 ~ 0.754` |
| 7 | branch step height at `x = 1/2` | `stair_diag` | definition-based, depends on the `prec` pair; see CI log |
| 8 | centre-point re-check, all 8 log segments | `stair_diag` | every segment reports `OK` |

Row 6 is measured against `mpmath` at 300 decimal digits and shows that the
tail factor `1/(1-x^2)` (evaluated as `4/3` at `x = 1/2`) leaves the bound
conservative by about 33%. Row 5 is an independent implementation of the
coefficient recursion; it is not an assertion on the same code path.

## Known limitations

A bound is only useful if its weakest link is visible, so these are stated
explicitly.

- **The direct branch's error is a magnitude estimate, not a rigorous bound.**
  It is the observed difference between two `prec` settings. Only the series
  branch carries an analytic bound.
- **There is a branch step at `x = 1/2`,** of the order of tens of orders of
  magnitude. It is an artifact of the `if x > SAFE` switch, not a property of
  `R(x)`, which is analytic across that point. The step height depends on the
  `prec` pair used for the direct branch (`stair_diag.py` uses `30` and `40`).
- **No comparison against Arb/FLINT has been performed.** This comparison is
  not done, and no claim of comparative accuracy is made.
- **`Phi(x) = R(x) - Si(x)` has infinitely many stationary points.** They are
  not enumerated here. `Phi(0+) = 0` and `Phi(inf) = -pi/2`, so no additive
  constant exists.
- **`R(x) < 0` is proved strictly only for `x >= 4`** (which follows from
  `x > 2*sqrt(3)`). For `x < 4` the sign is sampled, not proved.
- **Valid for `prec <= 30`.** With `K = 24` the series branch caps near
  `4.3e-56` at `x = 1/2`; raising `prec` past roughly 33 makes the series
  branch, not the direct branch, the limiting term.

## Reproducing these numbers

Every in-repository value in the table is produced by a command here.
`stair_diag.py` additionally writes its raw samples to `stair_diag.csv`, so the
figure can be regenerated rather than trusted.

The CI workflow runs all commands on GitHub's runners across three Python
versions and uploads the diagnostic outputs as build artifacts. The logs carry
the commit hash and a timestamp, so the results are reproducible from a machine
that is not the author's.

## License

MIT. See `LICENSE`.
