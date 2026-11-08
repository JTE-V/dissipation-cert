# Contributing

Thanks for looking at this repository. Bug reports and corrections are welcome.

## Reporting a problem

Open an issue and include:

- the command you ran, verbatim;
- the Python version (`python3 --version`);
- the full output, including the exit code;
- if the output is a number you believe is wrong, say what you expected instead.

If the problem is a wrong bound, a counterexample is the most useful thing you
can send. A single value of `x` where the printed bound is violated is worth
more than a description.

## Suggesting a change

This package aims to stay dependency-free. Pull requests that add a third-party
dependency will not be merged. Where an outside library is useful — for example
as an independent cross-check — the intent is to keep that check outside the
package, as a separate script, and to say so in the documentation.

## Tests

Run the test suite before opening a pull request:
