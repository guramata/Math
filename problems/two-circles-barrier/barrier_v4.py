"""
One-sided cut  B1(psi):

    Arc : circle 1,  phi in [-pi/2, pi/2 - psi]       (right semicircle minus a top cap)
    V   : x = -50,   h* <= y <= 100/cos(psi)
    H   : y = 0,     xa <= x <= xb

The design step sweeps every F-line missing the arc, records where it crosses
x = -50 and y = 0, and picks (h*, xa, xb) minimising |V| + |H| subject to covering
all of them.  Verification is by independent random sampling with a margin.
"""

import math
import random
import sys

R, A = 100.0, 50.0


def g(t):
    return R - A * abs(math.cos(t))


def arc_hit(theta, p, phi_lo, phi_hi, tol=1e-12):
    s = p + A * math.cos(theta)
    if abs(s) > R + tol:
        return False
    d = math.acos(max(-1.0, min(1.0, s / R)))
    for phi in (theta + d, theta - d):
        k = (phi - phi_lo) % (2 * math.pi) + phi_lo
        if phi_lo - tol <= k <= phi_hi + tol:
            return True
    return False


def seg_hit(theta, p, P, Q, tol=1e-9):
    c, s = math.cos(theta), math.sin(theta)
    f0 = P[0] * c + P[1] * s - p
    f1 = Q[0] * c + Q[1] * s - p
    return f0 * f1 <= 0.0 or abs(f0) <= tol or abs(f1) <= tol


def cross_x(theta, p, x0):
    s = math.sin(theta)
    return None if abs(s) < 1e-15 else (p - x0 * math.cos(theta)) / s


def cross_y(theta, p, y0):
    c = math.cos(theta)
    return None if abs(c) < 1e-15 else (p - y0 * math.sin(theta)) / c


def design(psi, n_t=2500, n_p=2500):
    lo, hi = -math.pi / 2, math.pi / 2 - psi
    top = R / math.cos(psi)
    rows = []
    for i in range(n_t + 1):
        t = math.pi * i / n_t
        gt = g(t)
        for j in range(n_p + 1):
            p = -gt + 2 * gt * j / n_p
            if arc_hit(t, p, lo, hi):
                continue
            yv = cross_x(t, p, -A)
            xh = cross_y(t, p, 0.0)
            rows.append((yv if yv is not None else -math.inf,
                         xh if xh is not None else math.nan))
    # scan h*: lines with yv in [h*, top] are covered by V; the rest must be covered by H
    rows.sort(key=lambda r: r[0])
    ys = [r[0] for r in rows]
    best = None
    cand = sorted(set([top] + [y for y in ys if 0 <= y <= top]))
    # prefix over rows sorted by yv ascending: lines with yv < h* -> need H
    # also lines with yv > top always need H
    always = [r[1] for r in rows if r[0] > top or math.isnan(r[1]) and False]
    if any(math.isnan(r[1]) and (r[0] < 0 or r[0] > top) for r in rows):
        return None  # a vertical line that V does not catch: design impossible
    tail_min = [math.inf] * (len(rows) + 1)
    tail_max = [-math.inf] * (len(rows) + 1)
    # H must contain xh of every row with yv > top  (compute once)
    hi_rows = [r[1] for r in rows if r[0] > top]
    base_lo = min(hi_rows) if hi_rows else math.inf
    base_hi = max(hi_rows) if hi_rows else -math.inf
    # prefix min/max of xh over rows with yv < h* (rows sorted by yv)
    pre_lo, pre_hi = [], []
    m1, m2 = math.inf, -math.inf
    for r in rows:
        if not math.isnan(r[1]):
            m1, m2 = min(m1, r[1]), max(m2, r[1])
        pre_lo.append(m1)
        pre_hi.append(m2)
    import bisect
    for hstar in cand:
        k = bisect.bisect_left(ys, hstar)          # rows[0..k-1] have yv < h*
        xa = min(base_lo, pre_lo[k - 1] if k > 0 else math.inf)
        xb = max(base_hi, pre_hi[k - 1] if k > 0 else -math.inf)
        lenH = 0.0 if xa == math.inf else (xb - xa)
        L = R * (hi - lo) + (top - hstar) + lenH
        if best is None or L < best[0]:
            best = (L, hstar, xa, xb)
    return best


def verify(psi, hstar, xa, xb, n_random=300_000, seed=3, margin=0.0):
    lo, hi = -math.pi / 2, math.pi / 2 - psi
    top = R / math.cos(psi)
    segs = [((-A, hstar + margin), (-A, top - margin))]
    if xa != math.inf:
        segs.append(((xa + margin, 0.0), (xb - margin, 0.0)))
    rng = random.Random(seed)
    bad = []
    for _ in range(n_random):
        t = rng.uniform(0.0, math.pi)
        p = rng.uniform(-g(t), g(t))
        if not (arc_hit(t, p, lo, hi) or any(seg_hit(t, p, P, Q) for P, Q in segs)):
            bad.append((t, p))
            if len(bad) >= 5:
                break
    return bad


if __name__ == "__main__":
    degs = [float(x) for x in sys.argv[1:]] or [5, 10, 15, 20, 25, 30, 35, 40, 45]
    print(" psi(deg)    length      h*        xa        xb   |V|     |H|    esc  esc(m=.05)")
    for deg in degs:
        psi = math.radians(deg)
        res = design(psi)
        if res is None:
            print(f"{deg:8.2f}   impossible")
            continue
        L, hstar, xa, xb = res
        top = R / math.cos(psi)
        b0 = verify(psi, hstar, xa, xb, n_random=150_000)
        b1 = verify(psi, hstar, xa, xb, n_random=150_000, margin=0.05)
        print(f"{deg:8.2f} {L:11.4f} {hstar:8.3f} {xa:9.3f} {xb:9.3f} {top-hstar:6.2f} {xb-xa:7.2f}   "
              f"{len(b0):2d}    {len(b1):2d}")
