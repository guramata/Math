"""Integer set-cover on the grid family (same setup as lp_barrier.py), HiGHS MILP with time limit."""
import math, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix, vstack

R, A = 100.0, 50.0
sp = float(sys.argv[1]); NT = int(sys.argv[2]); NP = int(sys.argv[3]); tlim = float(sys.argv[4]); out = sys.argv[5]
XMIN, XMAX, YMIN, YMAX = -160.0, 160.0, -115.0, 115.0
xs = np.arange(XMIN, XMAX + 1e-9, sp); ys = np.arange(YMIN, YMAX + 1e-9, sp)
nx, ny = len(xs), len(ys)
nodes = np.array([(x, y) for x in xs for y in ys])
dirs = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (1, 2), (2, -1), (1, -2)]
edges = []
for i in range(nx):
    for j in range(ny):
        for dx, dy in dirs:
            a, b = i + dx, j + dy
            if 0 <= a < nx and 0 <= b < ny:
                edges.append((i * ny + j, a * ny + b))
edges = np.array(edges); P = nodes[edges[:, 0]]; Q = nodes[edges[:, 1]]
lengths = np.linalg.norm(P - Q, axis=1)
ths = (np.arange(NT) + 0.5) * math.pi / NT
rows = []
for t in ths:
    gt = R - A * abs(math.cos(t))
    ps = -gt + (np.arange(NP) + 0.5) * 2 * gt / NP
    c, s = math.cos(t), math.sin(t)
    fP = P[:, 0] * c + P[:, 1] * s; fQ = Q[:, 0] * c + Q[:, 1] * s
    lo = np.minimum(fP, fQ); hi = np.maximum(fP, fQ)
    rows.append(csr_matrix(((lo[None, :] <= ps[:, None]) & (ps[:, None] <= hi[None, :])).astype(np.int8)))
Amat = vstack(rows).tocsr()
print(f"edges {len(edges)} lines {Amat.shape[0]}", flush=True)
t0 = time.time()
res = milp(lengths, constraints=LinearConstraint(Amat, lb=1, ub=np.inf), integrality=np.ones(len(edges)),
           bounds=Bounds(0, 1), options={"time_limit": tlim, "disp": True, "mip_rel_gap": 0.005})
print("status:", res.message, " time %.0fs" % (time.time() - t0))
x = np.round(res.x).astype(int) if res.x is not None else None
if x is None:
    sys.exit(1)
print(f"MILP length: {lengths[x == 1].sum():.3f}   (lower bound reported: {getattr(res, 'mip_dual_bound', float('nan')):.3f})")
np.save(out + ".npy", np.hstack([P[x == 1], Q[x == 1]]))
S = 2.2; W = (XMAX - XMIN) * S + 40; H = (YMAX - YMIN) * S + 40
px = lambda x_, y_: ((x_ - XMIN) * S + 20, (YMAX - y_) * S + 20)
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" viewBox="0 0 {W:.0f} {H:.0f}">',
       f'<rect width="{W:.0f}" height="{H:.0f}" fill="#fbfaf7"/>']
for cx in (-A, A):
    a = px(cx, 0)
    svg.append(f'<circle cx="{a[0]:.1f}" cy="{a[1]:.1f}" r="{R*S:.1f}" fill="none" stroke="#c3bdb0" stroke-dasharray="6 4"/>')
for k in np.where(x == 1)[0]:
    a, b = px(*P[k]), px(*Q[k])
    svg.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" stroke="#c2543a" stroke-width="3" stroke-linecap="round"/>')
svg.append(f'<text x="20" y="{H-8:.0f}" font-family="Georgia" font-size="14" fill="#22201d">MILP set cover, grid {sp:g}, length {lengths[x==1].sum():.1f} (sampled lines only)</text>')
svg.append('</svg>')
open(out, "w").write("\n".join(svg)); print("wrote", out)
