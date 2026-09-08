"""
Hybrid barrier  Y(psi):

    Arc1 : circle 1 (centre (-50,0)),  phi in [0, pi/2 - psi]          (upper-right quarter, pole end cut)
    Arc2 : circle 2 (centre ( 50,0)),  phi in [pi, 3pi/2 - psi]        (lower-left quarter, pole end cut)
    V+   : x = -50,   h* <= y <= 100/cos(psi)
    V-   : x = +50,  -100/cos(psi) <= y <= -h*        (180-degree rotation of V+)

Length  100(pi - 2 psi) + 2(100/cos psi - h*).
This script sweeps all F-lines missing both arcs, records where they cross the
two vertical lines, finds the smallest h* that catches all of them (or reports
lines that cross neither vertical line in the admissible height range), and
re-verifies by random sampling.
"""

import math
import random

R, A = 100.0, 50.0


def g(t):
    return R - A * abs(math.cos(t))


def arc_hit_general(theta, p, cx, cy, phi_lo, phi_hi, tol=1e-12):
    """Line x cos t + y sin t = p meets the arc centre (cx,cy), radius R, phi in [lo,hi]?"""
    s = p - (cx * math.cos(theta) + cy * math.sin(theta))
    if abs(s) > R + tol:
        return False
    d = math.acos(max(-1.0, min(1.0, s / R)))
    for phi in (theta + d, theta - d):
        # normalise phi into [phi_lo, phi_lo + 2pi)
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


def arcs_for(psi):
    return [(-A, 0.0, 0.0, math.pi / 2 - psi),
            (A, 0.0, math.pi, 3 * math.pi / 2 - psi)]


def design(psi, n_t=3000, n_p=3000):
    arcs = arcs_for(psi)
    top = R / math.cos(psi)
    worst = 0.0          # the largest h* we can afford = min over lines of the best crossing height
    hstar = top
    uncaught = []
    for i in range(n_t + 1):
        t = math.pi * i / n_t
        gt = g(t)
        for j in range(n_p + 1):
            p = -gt + 2 * gt * j / n_p
            if any(arc_hit_general(t, p, *a) for a in arcs):
                continue
            y1 = cross_x(t, p, -A)          # needs y1 in [h*, top]
            y2 = cross_x(t, p, A)           # needs y2 in [-top, -h*]
            cands = []
            if y1 is not None and 0 <= y1 <= top:
                cands.append(y1)
            if y2 is not None and -top <= y2 <= 0:
                cands.append(-y2)
            if not cands:
                uncaught.append((t, p, y1, y2))
                continue
            hstar = min(hstar, max(cands))
    return hstar, top, uncaught


def blocked(theta, p, arcs, segs):
    return any(arc_hit_general(theta, p, *a) for a in arcs) or \
        any(seg_hit(theta, p, P, Q) for P, Q in segs)


def verify(psi, hstar, n_random=300_000, seed=5, margin=0.0):
    arcs = arcs_for(psi)
    top = R / math.cos(psi)
    segs = [((-A, hstar - margin), (-A, top)), ((A, -hstar + margin), (A, -top))]
    rng = random.Random(seed)
    bad = []
    for _ in range(n_random):
        t = rng.uniform(0.0, math.pi)
        p = rng.uniform(-g(t), g(t))
        if not blocked(t, p, arcs, segs):
            bad.append((t, p))
            if len(bad) >= 5:
                break
    return bad


if __name__ == "__main__":
    print(" psi(deg)    length      h*     uncaught   escapes  escapes(margin .05)")
    for deg in [0, 5, 10, 15, 20, 25, 30, 35, 40]:
        psi = math.radians(deg)
        hstar, top, unc = design(psi, n_t=1500, n_p=1500)
        L = R * (math.pi - 2 * psi) + 2 * (top - hstar)
        bad0 = verify(psi, hstar, n_random=150_000)
        bad1 = verify(psi, hstar, n_random=150_000, margin=0.05)
        print(f"{deg:8.2f} {L:11.4f} {hstar:8.3f}   {len(unc):6d}   {len(bad0):3d}      {len(bad1):3d}   "
              f"{unc[:1]}")
