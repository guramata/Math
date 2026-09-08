"""
Set-cover LP heuristic for the two-circle barrier.

Candidate pieces: short segments joining grid points (8 directions).
Constraints: every sampled F-line must be hit by total weight >= 1.
Objective: minimise total length.  The LP relaxation value is a heuristic
lower bound for barriers built from these pieces (lines are sampled, so it is
not a rigorous bound); the fractional solution shows where length wants to be.

usage: python3 lp_barrier.py SPACING NTHETA NP [out.svg]
"""

import math
import sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack

R, A = 100.0, 50.0
sp = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
NT = int(sys.argv[2]) if len(sys.argv) > 2 else 90
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 100
out = sys.argv[4] if len(sys.argv) > 4 else None

XMIN, XMAX, YMIN, YMAX = -160.0, 160.0, -115.0, 115.0
xs = np.arange(XMIN, XMAX + 1e-9, sp)
ys = np.arange(YMIN, YMAX + 1e-9, sp)
nx, ny = len(xs), len(ys)
idx = lambda i, j: i * ny + j
nodes = np.array([(x, y) for x in xs for y in ys])

dirs = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (1, 2), (2, -1), (1, -2)]
edges = []
for i in range(nx):
    for j in range(ny):
        for dx, dy in dirs:
            a, b = i + dx, j + dy
            if 0 <= a < nx and 0 <= b < ny:
                edges.append((idx(i, j), idx(a, b)))
edges = np.array(edges)
P = nodes[edges[:, 0]]
Q = nodes[edges[:, 1]]
lengths = np.linalg.norm(P - Q, axis=1)
ne = len(edges)

# sampled F-lines: theta in (0,pi), p in (-g, g)
ths = (np.arange(NT) + 0.5) * math.pi / NT
rows = []
for t in ths:
    gt = R - A * abs(math.cos(t))
    ps = -gt + (np.arange(NP) + 0.5) * 2 * gt / NP
    c, s = math.cos(t), math.sin(t)
    fP = P[:, 0] * c + P[:, 1] * s          # (ne,)
    fQ = Q[:, 0] * c + Q[:, 1] * s
    # hit iff (fP - p)(fQ - p) <= 0
    lo = np.minimum(fP, fQ)
    hi = np.maximum(fP, fQ)
    M = (lo[None, :] <= ps[:, None]) & (ps[:, None] <= hi[None, :])
    rows.append(csr_matrix(M.astype(np.int8)))
Amat = vstack(rows).tocsr()
nl = Amat.shape[0]
print(f"grid {nx}x{ny}, edges {ne}, lines {nl}, nnz {Amat.nnz}", flush=True)

res = linprog(lengths, A_ub=-Amat, b_ub=-np.ones(nl), bounds=(0, 1),
              method="highs", options={"disp": False})
print("LP status:", res.message)
x = res.x
print(f"LP value (heuristic lower bound within grid family): {res.fun:.3f}")
sel = x > 1e-6
print(f"edges with weight>0: {sel.sum()},  weight>=0.5: {(x >= 0.5).sum()}, "
      f"integral length if rounding >=0.5: {lengths[x >= 0.5].sum():.2f}")

if out:
    S = 2.2
    W = (XMAX - XMIN) * S + 40
    H = (YMAX - YMIN) * S + 40
    px = lambda x, y: ((x - XMIN) * S + 20, (YMAX - y) * S + 20)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" viewBox="0 0 {W:.0f} {H:.0f}">',
           f'<rect width="{W:.0f}" height="{H:.0f}" fill="#fbfaf7"/>']
    for cx in (-A, A):
        a = px(cx, 0)
        svg.append(f'<circle cx="{a[0]:.1f}" cy="{a[1]:.1f}" r="{R*S:.1f}" fill="none" stroke="#c3bdb0" stroke-dasharray="6 4"/>')
    order = np.argsort(x)
    for k in order:
        if x[k] < 0.02:
            continue
        a, b = px(*P[k]), px(*Q[k])
        svg.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" '
                   f'stroke="#c2543a" stroke-width="{1+4*x[k]:.2f}" opacity="{min(1, 0.15+x[k]):.2f}"/>')
    svg.append(f'<text x="20" y="{H-8:.0f}" font-family="Georgia" font-size="14" fill="#22201d">'
               f'LP relaxation, grid {sp:g}, value {res.fun:.1f}; stroke weight = LP weight</text>')
    svg.append('</svg>')
    open(out, "w").write("\n".join(svg))
    print("wrote", out)
