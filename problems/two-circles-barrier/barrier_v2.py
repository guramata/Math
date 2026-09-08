"""
Improved barrier  B(psi)  for the two-circle problem (R = 100, centres +-50).

Start from the right semicircle S* of circle 1 (length 100 pi) and cut the two
arcs of angular size psi at its poles.  Replace them by three segments:

    V+ : vertical  x = -50,   h* <= y <= 100 / cos(psi)
    V- : mirror image of V+ in the x-axis
    H  : horizontal y = 0,   -50 - delta <= x <= -50 + 100 sin(psi)

h* and delta are the smallest values that catch every line of F missing the
retained arc; they are computed here by a dense sweep and then the whole
barrier is re-verified by independent random sampling.
"""

import math
import random

R, A = 100.0, 50.0


def g(t):
    return R - A * abs(math.cos(t))


# ------------------------------------------------------------- geometry
def seg_hit(theta, p, P, Q, tol=1e-9):
    c, s = math.cos(theta), math.sin(theta)
    f0 = P[0] * c + P[1] * s - p
    f1 = Q[0] * c + Q[1] * s - p
    return f0 * f1 <= 0.0 or abs(f0) <= tol or abs(f1) <= tol


def arc_hit(theta, p, phi_max, tol=1e-12):
    """Line meets the arc {C1 + R(cos phi, sin phi) : |phi| <= phi_max} ?"""
    s = p + A * math.cos(theta)
    if abs(s) > R + tol:
        return False
    d = math.acos(max(-1.0, min(1.0, s / R)))
    for phi in (theta + d, theta - d):
        phi = (phi + math.pi) % (2 * math.pi) - math.pi
        if abs(phi) <= phi_max + tol:
            return True
    return False


def cross_x(theta, p, x0):
    """y where the line meets the vertical line x = x0 (None if parallel)."""
    s = math.sin(theta)
    return None if abs(s) < 1e-15 else (p - x0 * math.cos(theta)) / s


def cross_y(theta, p, y0):
    c = math.cos(theta)
    return None if abs(c) < 1e-15 else (p - y0 * math.sin(theta)) / c


# ------------------------------------------------------------- design
def design(psi, n_t=3000, n_p=3000):
    """
    Sweep the F-lines that miss the retained arc and return the smallest
    (h*, delta) such that V+, V-, H catch all of them, minimising total length.
    """
    phi_max = math.pi / 2 - psi
    top = R / math.cos(psi)
    h_right = -A + R * math.sin(psi)
    need = []                      # lines caught by neither the V-top part nor H's [-50, h_right]
    for i in range(n_t + 1):
        t = math.pi * i / n_t
        gt = g(t)
        for j in range(n_p + 1):
            p = -gt + 2 * gt * j / n_p
            if arc_hit(t, p, phi_max):
                continue
            yv = cross_x(t, p, -A)
            xh = cross_y(t, p, 0.0)
            if xh is not None and -A <= xh <= h_right:
                continue
            if yv is None:                       # vertical line missing H: impossible in F
                need.append((0.0, -math.inf))
                continue
            need.append((abs(yv), xh if xh is not None else -math.inf))
    # every remaining line must satisfy |yv| >= h*  or  xh >= -50 - delta
    need.sort(key=lambda z: z[1])                # by x-crossing
    best = None
    # candidate deltas: each line's x-crossing defines a delta; h* = min |yv| over lines with smaller xh
    hs = [z[0] for z in need]
    # prefix minima of |yv| over lines sorted by xh ascending
    pref = []
    m = math.inf
    for h in hs:
        m = min(m, h)
        pref.append(m)
    for k in range(len(need) + 1):
        # lines 0..k-1 are covered by V (need h* <= pref[k-1]); lines k.. covered by H (delta = -50 - xh_k)
        hstar = pref[k - 1] if k > 0 else top
        delta = (-A - need[k][1]) if k < len(need) else 0.0
        if delta < 0:
            delta = 0.0
        if hstar > top:
            hstar = top
        L = R * (math.pi - 2 * psi) + 2 * (top - hstar) + (R * math.sin(psi) + delta)
        if best is None or L < best[0]:
            best = (L, hstar, delta)
    return best


def make_barrier(psi, hstar, delta):
    top = R / math.cos(psi)
    segs = [
        ((-A, hstar), (-A, top)),
        ((-A, -hstar), (-A, -top)),
        ((-A - delta, 0.0), (-A + R * math.sin(psi), 0.0)),
    ]
    return segs, math.pi / 2 - psi


def blocked(theta, p, segs, phi_max):
    return arc_hit(theta, p, phi_max) or any(seg_hit(theta, p, P, Q) for P, Q in segs)


def verify(psi, hstar, delta, n_random=300_000, seed=11, margin=0.0):
    segs, phi_max = make_barrier(psi, hstar - margin, delta + margin)
    rng = random.Random(seed)
    bad = []
    for _ in range(n_random):
        t = rng.uniform(0.0, math.pi)
        p = rng.uniform(-g(t), g(t))
        if not blocked(t, p, segs, phi_max):
            bad.append((t, p))
            if len(bad) >= 5:
                break
    return bad


if __name__ == "__main__":
    print(" psi(deg)    length      h*      delta   escapes(no margin)  escapes(margin 0.05)")
    for deg in [4, 6, 8, 10, 11, 12, 13, 14, 16, 18, 20]:
        psi = math.radians(deg)
        L, hstar, delta = design(psi)
        bad0 = verify(psi, hstar, delta)
        bad1 = verify(psi, hstar, delta, margin=0.05)
        print(f"{deg:8.2f} {L:11.4f} {hstar:8.3f} {delta:8.3f}   {len(bad0):3d}                {len(bad1):3d}")
