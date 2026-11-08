#!/usr/bin/env python3
"""test_all.py — assert the values reported by CI, using only the standard library.

Runs the two scripts as subprocesses and checks that their printed output
contains the values the README documents. Deliberately avoids importing
internal functions, so it stays valid across refactors.

Usage: python3 test_all.py
Exit code 0 on success, 1 on any failed check.
"""
import subprocess
import sys

FAILS = []


def run(script, *args):
    p = subprocess.run([sys.executable, script, *args],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def check(label, cond, detail=""):
    print(("  ok    " if cond else "  FAIL  ") + label + ("   " + detail if detail else ""))
    if not cond:
        FAILS.append(label)


print("exact")
rc, out = run("dissipation_cert.py", "exact")
check("exit code 0", rc == 0)
for c in ("-6", "-3/10", "-11/1400", "-17/126000"):
    check("coefficient %s reported" % c, c in out)

print("selfcheck")
rc, out = run("dissipation_cert.py", "selfcheck")
check("exit code 0", rc == 0)
check("first dropped term at x=1/2", "3.19e-56" in out)
check("two-path collision at x=1/2", "1.616e-33" in out)
check("Si(pi) to 30 digits",
      "1.851937051982466170361053370157" in out)

print("stair_diag")
rc, out = run("stair_diag.py")
check("exit code 0", rc == 0)
check("|H_24| = 1.123705e-42", "1.123705e-42" in out)
check("series bound at x=1/2, tail included", "4.258e-56" in out)
check("series bound at x=1/2, dropped term only", "3.194e-56" in out)
check("all 8 segment centres below the noise floor",
      out.count("已在噪声") >= 8 or out.count("OK") >= 8)

print()
if FAILS:
    print("%d check(s) failed:" % len(FAILS))
    for f in FAILS:
        print("  - " + f)
    sys.exit(1)
print("all checks passed")
