"""
Two-circle barrier problem  (R = 100, centre distance d = 100).

C1 = (-50, 0), C2 = (50, 0), both of radius R = 100.

F  = family of lines meeting BOTH disks.
Line  L(theta, p):  x cos(theta) + y sin(theta) = p.

Signed distances of the centres from the line:
    s = p + 50 cos(theta)      (from C1)
    t = p - 50 cos(theta)      (from C2)

L in F  <=>  |s| <= 100 and |t| <= 100  <=>  |p| <= 100 - 50|cos theta| =: g(theta).

This script checks, numerically:
  1. the Cauchy-Crofton measure  mu(F) = 200(pi - 1);
  2. the lower bound  inf length >= mu(F)/2 = 100(pi - 1);
  3. that the right semicircle of C1 blocks every line of F  (upper bound 100 pi);
  4. the measure of F-lines that miss the lens D1 ^ D2  ( = 200pi/3 - 200 );
  5. that the two tip segments block exactly those lens-missing lines.
No third-party packages required.
"""

import math
import random

R = 100.0
A = 50.0            # half distance between centres
C1 = (-A, 0.0)
C2 = (A, 0.0)


def g(theta):
    """Half-width of F in the direction theta."""
    return R - A * abs(math.cos(theta))


def in_F(theta, p):
    s = p + A * math.cos(theta)
    t = p - A * math.cos(theta)
    return abs(s) <= R and abs(t) <= R


# ---------------------------------------------------------------- 1 & 2
def measure_F(n=2_000_000):
    """mu(F) = int_0^pi 2 g(theta) d(theta), by the trapezoid rule."""
    tot = 0.0
    for i in range(n + 1):
        th = math.pi * i / n
        w = 0.5 if i in (0, n) else 1.0
        tot += w * 2.0 * g(th)
    return tot * math.pi / n


# ---------------------------------------------------------------- 3
def hits_right_semicircle(theta, p, tol=1e-12):
    """
    Right semicircle of C1:  P(phi) = C1 + R(cos phi, sin phi),  phi in [-pi/2, pi/2].
    Substituting into the line equation:
        -A cos(theta) + R cos(phi - theta) = p
    so a hit exists iff  cos(phi - theta) = s/R  for some phi in [-pi/2, pi/2].
    """
    s = p + A * math.cos(theta)
    if abs(s) > R + tol:
        return False
    c = max(-1.0, min(1.0, s / R))
    d = math.acos(c)
    for phi in (theta + d, theta - d):
        # reduce to (-pi, pi]
        phi = (phi + math.pi) % (2 * math.pi) - math.pi
        if -math.pi / 2 - tol <= phi <= math.pi / 2 + tol:
            return True
    return False


def check_semicircle(trials=2_000_000, seed=1):
    rng = random.Random(seed)
    bad = []
    for _ in range(trials):
        th = rng.uniform(0.0, math.pi)
        p = rng.uniform(-g(th), g(th))
        if not hits_right_semicircle(th, p):
            bad.append((th, p))
            if len(bad) > 5:
                break
    # deterministic sweep as well
    N, M = 2000, 2000
    for i in range(N + 1):
        th = math.pi * i / N
        for j in range(M + 1):
            p = -g(th) + 2 * g(th) * j / M
            if not hits_right_semicircle(th, p):
                bad.append((th, p))
                if len(bad) > 5:
                    return bad
    return bad


# ---------------------------------------------------------------- 4 & 5
LENS_TIP = A * math.sqrt(3.0)          # 50*sqrt(3) = 86.6025...


def h_lens(theta):
    """Support function of the lens D1 ^ D2 (symmetric in theta -> theta+pi)."""
    c = abs(math.cos(theta))
    if c >= 0.5:                       # support point on a circular arc
        return R - A * c
    return LENS_TIP * abs(math.sin(theta))


def measure_F_missing_lens(n=2_000_000):
    """mu({ lines in F that miss the lens }) = 2 * int_0^pi (g - h_lens)_+ dtheta."""
    tot = 0.0
    for i in range(n + 1):
        th = math.pi * i / n
        w = 0.5 if i in (0, n) else 1.0
        tot += w * 2.0 * max(0.0, g(th) - h_lens(th))
    return tot * math.pi / n


def tip_segments_block_lens_missers(N=4000, M=4000):
    """
    Claim: every F-line missing the lens crosses the y-axis at height
    LENS_TIP < |y| <= R, i.e. it meets one of the two segments
        {0} x [ LENS_TIP, R ]   and   {0} x [ -R, -LENS_TIP ].
    """
    bad = []
    for i in range(1, N):
        th = math.pi * i / N
        lo, hi = h_lens(th), g(th)
        if hi <= lo + 1e-12:
            continue
        st = math.sin(th)
        if abs(st) < 1e-12:
            continue
        for j in range(M + 1):
            p = lo + (hi - lo) * j / M
            y = p / st                      # crossing height of the y-axis
            if not (LENS_TIP - 1e-9 <= abs(y) <= R + 1e-9):
                bad.append((th, p, y))
                if len(bad) > 5:
                    return bad
    return bad


if __name__ == "__main__":
    muF = measure_F()
    print("mu(F)              = %.9f   (exact 200(pi-1) = %.9f)"
          % (muF, 200 * (math.pi - 1)))
    print("Crofton lower bound= %.9f   (exact 100(pi-1) = %.9f)"
          % (muF / 2, 100 * (math.pi - 1)))

    bad = check_semicircle()
    print("right semicircle: unblocked F-lines found:", len(bad))
    print("upper bound        = %.9f   (exact 100 pi)" % (math.pi * R))

    muM = measure_F_missing_lens()
    print("mu(F \\ lines meeting lens) = %.9f   (exact 200pi/3 - 200 = %.9f)"
          % (muM, 200 * math.pi / 3 - 200))

    bad2 = tip_segments_block_lens_missers()
    print("tip segments: unblocked lens-missing lines:", len(bad2))
    print("tip segment total length = %.9f   (exact 200 - 100 sqrt 3)"
          % (2 * (R - LENS_TIP)))
