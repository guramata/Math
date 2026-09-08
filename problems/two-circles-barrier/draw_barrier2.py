"""
Draws the improved barrier B(psi): the right semicircle of circle 1 with both
pole caps cut, plus three segments.  Writes barrier2.svg.  Standard library only.
"""

import math

R, A = 100.0, 50.0
C1 = (-A, 0.0)
PSI_DEG = 9.0            # cap half-angle
H_STAR = 98.4            # bottom of the vertical segments (conservative)
DELTA = 1.6              # left extension of the axis segment (conservative)

psi = math.radians(PSI_DEG)
TOP = R / math.cos(psi)
XR = -A + R * math.sin(psi)
LEN = R * (math.pi - 2 * psi) + 2 * (TOP - H_STAR) + (R * math.sin(psi) + DELTA)

XMIN, XMAX, YMIN, YMAX = -170.0, 170.0, -118.0, 118.0
SCALE = 2.35
PAD_L = PAD_T = 26.0
PAD_B = 118.0
W = (XMAX - XMIN) * SCALE + 2 * PAD_L
H = (YMAX - YMIN) * SCALE + PAD_T + PAD_B
INK, MUTED, HAIR, PAPER, LENS = "#22201d", "#8a8478", "#c3bdb0", "#fbfaf7", "#ecebe4"
RAY, BAR, CUT = "#8fa9c9", "#c2543a", "#e5b8ad"


def px(x, y):
    return ((x - XMIN) * SCALE + PAD_L, (YMAX - y) * SCALE + PAD_T)


def clip(theta, p):
    ux, uy = math.cos(theta), math.sin(theta)
    ox, oy = p * ux, p * uy
    dx, dy = -uy, ux
    lo, hi = -1e9, 1e9
    for d, o, a, b in ((dx, ox, XMIN, XMAX), (dy, oy, YMIN, YMAX)):
        if abs(d) < 1e-12:
            if not (a <= o <= b):
                return None
            continue
        t0, t1 = (a - o) / d, (b - o) / d
        lo, hi = max(lo, min(t0, t1)), min(hi, max(t0, t1))
    return None if lo >= hi else ((ox + lo * dx, oy + lo * dy), (ox + hi * dx, oy + hi * dy))


def line_from_points(P, Q):
    dx, dy = Q[0] - P[0], Q[1] - P[1]
    n = math.hypot(dx, dy)
    ux, uy = -dy / n, dx / n
    theta = math.atan2(uy, ux)
    p = P[0] * ux + P[1] * uy
    if theta < 0:
        theta += math.pi; p = -p
    return theta, p


# lines that used to be caught only by the cut caps, one per replacement segment
pole = (-A, R)
p0 = (-A + R * math.sin(psi), R * math.cos(psi))
RAYS = [
    line_from_points((-A + 3, R - 0.05), (-A + 14, R * math.cos(psi) + 0.6)),   # cap chord -> V+ above y=100
    line_from_points((-A, 99.3), (-A - 90, 60)),                                 # near-horizontal cap-left line -> V+
    line_from_points((-A + 9, 99.6), (-A + 6, -99.8)),                           # cap-to-cap near-vertical -> H
    line_from_points((-A - 0.6, 0), (-A + 8, R * math.cos(psi) + 1.2)),          # just left of C1 -> H (delta part)
    (math.pi / 2, 40.0),
    (0.55, 20.0),
]

o = []
w = o.append
w(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" viewBox="0 0 {W:.0f} {H:.0f}" font-family="Georgia,serif">')
w(f'<rect width="{W:.0f}" height="{H:.0f}" fill="{PAPER}"/>')
w(f'<clipPath id="d1"><circle cx="{px(*C1)[0]:.2f}" cy="{px(*C1)[1]:.2f}" r="{R*SCALE:.2f}"/></clipPath>')
w(f'<circle cx="{px(A,0)[0]:.2f}" cy="{px(A,0)[1]:.2f}" r="{R*SCALE:.2f}" fill="{LENS}" clip-path="url(#d1)"/>')
for c in (C1, (A, 0.0)):
    cx, cy = px(*c)
    w(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{R*SCALE:.2f}" fill="none" stroke="{HAIR}" stroke-width="1.6" stroke-dasharray="7 5"/>')

for th, p in RAYS:
    seg = clip(th, p)
    if seg:
        a, b = px(*seg[0]), px(*seg[1])
        w(f'<line x1="{a[0]:.2f}" y1="{a[1]:.2f}" x2="{b[0]:.2f}" y2="{b[1]:.2f}" stroke="{RAY}" stroke-width="1.5" opacity="0.85"/>')

# cut caps, drawn faint
for sgn in (1, -1):
    a = px(-A, sgn * R)
    b = px(*p0) if sgn > 0 else px(p0[0], -p0[1])
    sweep = 1 if sgn > 0 else 0
    w(f'<path d="M {a[0]:.2f} {a[1]:.2f} A {R*SCALE:.2f} {R*SCALE:.2f} 0 0 {sweep} {b[0]:.2f} {b[1]:.2f}" fill="none" stroke="{CUT}" stroke-width="3" stroke-dasharray="4 5"/>')

# retained arc
a, b = px(*p0), px(p0[0], -p0[1])
w(f'<path d="M {a[0]:.2f} {a[1]:.2f} A {R*SCALE:.2f} {R*SCALE:.2f} 0 0 1 {b[0]:.2f} {b[1]:.2f}" fill="none" stroke="{BAR}" stroke-width="5.5" stroke-linecap="round"/>')
# segments
for P, Q in (((-A, H_STAR), (-A, TOP)), ((-A, -H_STAR), (-A, -TOP)), ((-A - DELTA, 0.0), (XR, 0.0))):
    a, b = px(*P), px(*Q)
    w(f'<line x1="{a[0]:.2f}" y1="{a[1]:.2f}" x2="{b[0]:.2f}" y2="{b[1]:.2f}" stroke="{BAR}" stroke-width="5.5" stroke-linecap="round"/>')

for c, lab in ((C1, "C₁"), ((A, 0.0), "C₂")):
    cx, cy = px(*c)
    w(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="3.2" fill="{INK}"/>')
    w(f'<text x="{cx:.2f}" y="{cy+20:.2f}" fill="{INK}" font-size="15" text-anchor="middle" font-style="italic">{lab}</text>')
for P, lab, dx, dy in (((-A, TOP), "V₊", -16, 4), ((-A, -TOP), "V₋", -16, 8), ((XR, 0), "H", 12, 5)):
    a = px(*P)
    w(f'<text x="{a[0]+dx:.2f}" y="{a[1]+dy:.2f}" fill="{BAR}" font-size="14" font-style="italic">{lab}</text>')

y0 = H - PAD_B + 26
w(f'<text x="{PAD_L:.0f}" y="{y0:.0f}" fill="{INK}" font-size="17">Improved barrier: the right half of the left circle with caps of {PSI_DEG:g}° cut at both poles, plus three segments.</text>')
w(f'<text x="{PAD_L:.0f}" y="{y0+24:.0f}" fill="{MUTED}" font-size="15">V±: x = −50, {H_STAR:g} ≤ |y| ≤ 100/cos ψ = {TOP:.2f}.   H: y = 0, −50 − {DELTA:g} ≤ x ≤ −50 + 100 sin ψ = {XR:.2f}.</text>')
w(f'<text x="{PAD_L:.0f}" y="{y0+48:.0f}" fill="{MUTED}" font-size="15">Total length {LEN:.2f}, versus 100π = 314.16 for the full semicircle. Thin lines: samples that the cut caps used to catch.</text>')
w(f'<line x1="{PAD_L:.0f}" y1="{y0+68:.0f}" x2="{PAD_L+52:.0f}" y2="{y0+68:.0f}" stroke="{BAR}" stroke-width="5.5" stroke-linecap="round"/>')
w(f'<text x="{PAD_L+62:.0f}" y="{y0+73:.0f}" fill="{MUTED}" font-size="14">barrier</text>')
w(f'<line x1="{PAD_L+140:.0f}" y1="{y0+68:.0f}" x2="{PAD_L+192:.0f}" y2="{y0+68:.0f}" stroke="{CUT}" stroke-width="3" stroke-dasharray="4 5"/>')
w(f'<text x="{PAD_L+202:.0f}" y="{y0+73:.0f}" fill="{MUTED}" font-size="14">removed caps</text>')
w('</svg>')
open(__file__.rsplit("/", 1)[0] + "/barrier2.svg", "w").write("\n".join(o))
print(f"wrote barrier2.svg  length={LEN:.3f}")
