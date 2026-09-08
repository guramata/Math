"""
Draws the shortest barrier known for the two-circle problem:
the half of the left circle lying in x >= -50, of length 100*pi.

Writes barrier.svg next to this script.  Pure standard library.
"""

import math

R, A = 100.0, 50.0
C1, C2 = (-A, 0.0), (A, 0.0)

# ---------------------------------------------------------------- canvas
XMIN, XMAX, YMIN, YMAX = -170.0, 170.0, -118.0, 118.0
SCALE = 2.35
PAD_L = PAD_T = 26.0
PAD_B = 96.0
W = (XMAX - XMIN) * SCALE + 2 * PAD_L
H = (YMAX - YMIN) * SCALE + PAD_T + PAD_B

INK      = "#22201d"
MUTED    = "#8a8478"
HAIR     = "#c3bdb0"
PAPER    = "#fbfaf7"
LENS     = "#ecebe4"
RAY      = "#8fa9c9"
BARRIER  = "#c2543a"


def px(x, y):
    return ((x - XMIN) * SCALE + PAD_L, (YMAX - y) * SCALE + PAD_T)


def clip(theta, p):
    """Line x cos t + y sin t = p, clipped to the drawing box."""
    ux, uy = math.cos(theta), math.sin(theta)
    ox, oy = p * ux, p * uy            # foot of the perpendicular
    dx, dy = -uy, ux                   # direction along the line
    lo, hi = -1e9, 1e9
    for d, o, a, b in ((dx, ox, XMIN, XMAX), (dy, oy, YMIN, YMAX)):
        if abs(d) < 1e-12:
            if not (a <= o <= b):
                return None
            continue
        t0, t1 = (a - o) / d, (b - o) / d
        lo, hi = max(lo, min(t0, t1)), min(hi, max(t0, t1))
    if lo >= hi:
        return None
    return ((ox + lo * dx, oy + lo * dy), (ox + hi * dx, oy + hi * dy))


def arc_hits(theta, p):
    """Where the line meets the barrier arc phi in [-pi/2, pi/2]."""
    s = p + A * math.cos(theta)                 # signed distance from C1
    if abs(s) > R:
        return []
    d = math.acos(max(-1.0, min(1.0, s / R)))
    out = []
    for phi in (theta + d, theta - d):
        phi = (phi + math.pi) % (2 * math.pi) - math.pi
        if -math.pi / 2 <= phi <= math.pi / 2:
            out.append((C1[0] + R * math.cos(phi), C1[1] + R * math.sin(phi)))
    return out


# lines meeting both circles, chosen to show the variety of incidences
RAYS = [
    (math.pi / 2,  99.0),      # nearly the top common tangent
    (math.pi / 2,  55.0),
    (math.pi / 2, -82.0),
    (0.0,        -49.0),       # nearly the vertical tangent x = -50
    (0.0,         44.0),
    (0.55,        20.0),
    (2.35,       -30.0),
    (1.15,       -70.0),
]

out = []
w = out.append
w(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" '
  f'viewBox="0 0 {W:.0f} {H:.0f}" font-family="Georgia,serif">')
w(f'<rect width="{W:.0f}" height="{H:.0f}" fill="{PAPER}"/>')

# lens
w(f'<clipPath id="d1"><circle cx="{px(*C1)[0]:.2f}" cy="{px(*C1)[1]:.2f}" '
  f'r="{R*SCALE:.2f}"/></clipPath>')
w(f'<circle cx="{px(*C2)[0]:.2f}" cy="{px(*C2)[1]:.2f}" r="{R*SCALE:.2f}" '
  f'fill="{LENS}" clip-path="url(#d1)"/>')

# the two circles
for c in (C1, C2):
    cx, cy = px(*c)
    w(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{R*SCALE:.2f}" fill="none" '
      f'stroke="{HAIR}" stroke-width="1.6" stroke-dasharray="7 5"/>')

# sample lines and their incidences with the barrier
dots = []
for th, p in RAYS:
    seg = clip(th, p)
    if seg is None:
        continue
    (x0, y0), (x1, y1) = seg
    a, b = px(x0, y0), px(x1, y1)
    w(f'<line x1="{a[0]:.2f}" y1="{a[1]:.2f}" x2="{b[0]:.2f}" y2="{b[1]:.2f}" '
      f'stroke="{RAY}" stroke-width="1.5" opacity="0.85"/>')
    dots.extend(arc_hits(th, p))

# the barrier: right half of circle 1
top, bot = px(C1[0], R), px(C1[0], -R)
w(f'<path d="M {top[0]:.2f} {top[1]:.2f} A {R*SCALE:.2f} {R*SCALE:.2f} 0 0 1 '
  f'{bot[0]:.2f} {bot[1]:.2f}" fill="none" stroke="{BARRIER}" '
  f'stroke-width="5.5" stroke-linecap="round"/>')

for (x, y) in dots:
    cx, cy = px(x, y)
    w(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="3.6" fill="{BARRIER}"/>')

# centres, and the 100 m between them
for c, lab in ((C1, "C₁"), (C2, "C₂")):
    cx, cy = px(*c)
    w(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="3.2" fill="{INK}"/>')
    w(f'<text x="{cx:.2f}" y="{cy+20:.2f}" fill="{INK}" font-size="15" '
      f'text-anchor="middle" font-style="italic">{lab}</text>')
a, b = px(*C1), px(*C2)
w(f'<line x1="{a[0]:.2f}" y1="{a[1]:.2f}" x2="{b[0]:.2f}" y2="{b[1]:.2f}" '
  f'stroke="{INK}" stroke-width="1.3"/>')
w(f'<text x="{(a[0]+b[0])/2:.2f}" y="{a[1]-9:.2f}" fill="{INK}" font-size="14" '
  f'text-anchor="middle">100</text>')

# poles of the barrier
for y, dy, lab in ((R, -12, "(\u221250, 100)"), (-R, 22, "(\u221250, \u2212100)")):
    cx, cy = px(C1[0], y)
    w(f'<text x="{cx:.2f}" y="{cy+dy:.2f}" fill="{BARRIER}" font-size="13" '
      f'text-anchor="middle">{lab}</text>')

# caption
y0 = H - PAD_B + 26
w(f'<text x="{PAD_L:.0f}" y="{y0:.0f}" fill="{INK}" font-size="17">'
  f'Shortest barrier known: the half of the left circle with x ≥ −50.'
  f'</text>')
w(f'<text x="{PAD_L:.0f}" y="{y0+26:.0f}" fill="{MUTED}" font-size="15">'
  f'Length 100π = 314.159… Every line meeting both circles crosses it. '
  f'The Crofton bound 100(π−1) = 214.159… is a strict lower bound.</text>')
w(f'<line x1="{PAD_L:.0f}" y1="{y0+44:.0f}" x2="{PAD_L+52:.0f}" y2="{y0+44:.0f}" '
  f'stroke="{BARRIER}" stroke-width="5.5" stroke-linecap="round"/>')
w(f'<text x="{PAD_L+62:.0f}" y="{y0+49:.0f}" fill="{MUTED}" font-size="14">barrier</text>')
w(f'<line x1="{PAD_L+140:.0f}" y1="{y0+44:.0f}" x2="{PAD_L+192:.0f}" y2="{y0+44:.0f}" '
  f'stroke="{RAY}" stroke-width="1.5"/>')
w(f'<text x="{PAD_L+202:.0f}" y="{y0+49:.0f}" fill="{MUTED}" font-size="14">'
  f'sample lines meeting both circles</text>')
w('</svg>')

path = __file__.rsplit("/", 1)[0] + "/barrier.svg"
open(path, "w").write("\n".join(out))
print("wrote", path, "-", len(dots), "incidences marked")
