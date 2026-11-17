"""GIDEON round 1: six artists, six marks. 64x64 viewBox. Family: AMR ground + seal."""
import math, random
import numpy as np
from brush import stroke, dot, poly, arc, Layer

AMR_OPS = (0.186, 0.384, 0.62)


def amr(L, r1, r2, seed=7):
    """Three-level brushed AMR quadtree. r1 refined at h=8 (aligned to 16), r2 at h=4 (aligned to 8)."""
    rng = random.Random(seed)
    def line(p, q, op, w):
        L.add("mesh", stroke([p, q], w, seed=rng.randint(0, 10**6), step=0.9, rise=0.05, tail=0.1,
                             tail_min=0.55, press=0.1, wobble=0.18), op)
    for v in (16, 32, 48):
        line((0.3, v + rng.uniform(-.2, .2)), (63.7, v + rng.uniform(-.2, .2)), AMR_OPS[0], 0.42)
        line((v + rng.uniform(-.2, .2), 0.3), (v + rng.uniform(-.2, .2), 63.7), AMR_OPS[0], 0.42)
    (a0, b0, a1, b1) = r1
    for v in range(a0 + 8, a1, 16):
        line((v, b0), (v, b1), AMR_OPS[1], 0.36)
    for v in range(b0 + 8, b1, 16):
        line((a0, v), (a1, v), AMR_OPS[1], 0.36)
    (a0, b0, a1, b1) = r2
    for v in range(a0 + 4, a1, 8):
        line((v, b0), (v, b1), AMR_OPS[2], 0.3)
    for v in range(b0 + 4, b1, 8):
        line((a0, v), (a1, v), AMR_OPS[2], 0.3)


def seal(L, x=52.7, y=53.6, s=9.3, seed=3):
    """Carved seal, lower right: a square in the accent with a rectilinear G cut out (evenodd)."""
    rng = random.Random(seed)
    j = lambda: rng.uniform(-0.18, 0.18)
    outer = [(x + j(), y + j()), (x + s + j(), y - 0.12 + j()), (x + s + 0.1 + j(), y + s + j()), (x + j() - 0.05, y + s + 0.1 + j())]
    g = [(0.2, 0.2), (0.8, 0.2), (0.8, 0.33), (0.33, 0.33), (0.33, 0.67), (0.67, 0.67), (0.67, 0.6),
         (0.5, 0.6), (0.5, 0.48), (0.8, 0.48), (0.8, 0.8), (0.2, 0.8)]
    inner = [(x + u * s + rng.uniform(-.08, .08), y + v * s + rng.uniform(-.08, .08)) for u, v in g]
    L.add("acc", [np.array(outer), np.array(inner)], None, evenodd=True)


# ---------------------------------------------------------------- A · Karesansui
def mark_A():
    """The garden monk: five places on raked gravel; the measured one holds only ripples.
    Stones sized and toned by the 4th-order first derivative weights (+1, -8, 0, +8, -1)/12."""
    L = Layer()
    amr(L, (16, 16, 48, 48), (24, 24, 40, 40))
    xs = [8.5, 20.0, 31.5, 43.0, 54.5]
    yS = 31.0
    radii = [1.7, 3.5, 0.0, 3.5, 1.7]
    signs = [+1, -1, 0, +1, -1]
    ring_r = [r + 2.4 for r in radii]
    ring_r[2] = 2.6

    def bump(x, y0, side):
        y = y0
        for xs_, rr in zip(xs, ring_r):
            reach = rr + 3.6
            if abs(y0 - yS) < reach + 2:
                y += side * max(0, (reach + 1.0) - abs(y0 - yS)) * math.exp(-((x - xs_) / (rr + 2.2)) ** 2)
        return y
    # rake lines, three tines each
    sd = 20
    for y0, side in ((15.5, -1), (22.5, -1), (39.5, +1), (46.5, +1)):
        for tine in (-0.75, 0, 0.75):
            pts = [(x, bump(x, y0 + tine, side)) for x in np.linspace(3.5, 60.5, 40)]
            if y0 > 40:
                pts = [(x, y) for x, y in pts if x < 49.5]  # leave room for the seal
            sd += 1
            L.add("ink", stroke(pts, 0.42, seed=sd, step=0.4, rise=0.04, tail=0.08, tail_min=0.4, wobble=0.12))
    # ripples around each place
    for i, (x, r, rr, sg) in enumerate(zip(xs, radii, ring_r, signs)):
        role = "acc" if i == 2 else "ink"
        for k, dr in enumerate((0, 1.5, 3.0)):
            if i in (0, 4) and k == 2:
                continue
            a0 = random.Random(i * 10 + k).uniform(0, 6.28)
            pts = arc(x, yS, rr + dr, a0, a0 + 2 * math.pi * 0.96, 48)
            L.add(role, stroke(pts, 0.5 if role == "ink" else 0.62, seed=100 + i * 10 + k, step=0.3, rise=0.05, tail=0.1, tail_min=0.35, wobble=0.1))
        if r > 0:
            col = "ink" if sg > 0 else "mesh"
            L.add(col, dot(x, yS, r, seed=40 + i, squash=0.78, angle=0.3 * sg, rough=0.12), None if col == "ink" else 0.5)
    seal(L)
    return L


# ---------------------------------------------------------------- B · Take no fushi
def mark_B():
    """The bamboo painter: a stalk grows in equal intervals; the faint second stalk has its
    nodes halfway (the faces). The sprig grows from the node being measured."""
    L = Layer()
    amr(L, (32, 16, 64, 48), (32, 24, 48, 40))
    # faint stalk behind, nodes staggered by h/2
    def stalk(x0, y0, x1, y1, nseg, w, role, op, seed, nodes_off=0.0, bend=1.6):
        P0, P1 = np.array([x0, y0]), np.array([x1, y1])
        def at(t):
            p = P0 + (P1 - P0) * t
            return (p[0] + bend * math.sin(math.pi * t), p[1])
        cuts = [(k + nodes_off) / nseg for k in range(nseg + 1)]
        cuts = [min(max(c, 0), 1) for c in cuts]
        nodes = []
        for k in range(len(cuts) - 1):
            a, b = cuts[k] + 0.012, cuts[k + 1] - 0.012
            if b - a < 0.03:
                continue
            seg = [at(a + (b - a) * u) for u in np.linspace(0, 1, 6)]
            bp = lambda t: (0.82 + 0.18 * np.clip(t / 0.12, 0, 1)) * (1 + 0.1 * np.exp(-((t - 0.04) / 0.05) ** 2) + 0.12 * np.exp(-((t - 0.97) / 0.04) ** 2))
            L.add(role, stroke(seg, w, seed=seed + k, step=0.3, profile=bp, wobble=0.04, bristles=5, kasure=0.0), op)
            nodes.append(cuts[k + 1])
        return at, nodes
    at_f, nodes_f = stalk(20.5, 60, 23.5, 6, 5, 2.4, "mesh", 0.5, 300, nodes_off=0.5, bend=-1.2)
    for t in nodes_f[:-1]:
        x, y = at_f(t)
        L.add("mesh", stroke([(x - 2.2, y + 0.3), (x, y - 0.25), (x + 2.2, y + 0.3)], 0.8, seed=int(t * 1000), step=0.25, tail_min=0.3), 0.5)
    at, nodes = stalk(31.0, 61, 35.5, 3, 5, 4.4, "ink", None, 200)
    for i, t in enumerate(nodes[:-1]):
        x, y = at(t)
        L.add("ink", stroke([(x - 3.6, y + 0.6), (x - 1, y - 0.3), (x + 1.4, y - 0.4), (x + 3.5, y + 0.5)], 1.15,
                            seed=500 + i, step=0.25, rise=0.15, tail=0.3, tail_min=0.25, press=0.3))
    # sprig at the middle node (index 2 of 4 interior nodes -> the measured one)
    nx, ny = at(nodes[2])
    leaf = lambda t: np.sin(np.pi * np.clip(t, 0, 1)) ** 0.75 * (1 - 0.15 * t) + 0.04
    for k, (dx, dy, w) in enumerate(((16, -7, 3.6), (17, 2.5, 3.2), (11, -13, 2.8))):
        cx, cy = nx + 2.5, ny - 0.5
        pts = [(cx, cy), (cx + dx * 0.45, cy + dy * 0.35 - 1.2), (cx + dx, cy + dy)]
        L.add("acc", stroke(pts, w, seed=700 + k, step=0.25, profile=leaf, wobble=0.05))
    L.add("acc", stroke([(nx + 1.5, ny - 0.2), (nx + 4.0, ny - 1.0)], 0.9, seed=790, step=0.2, tail_min=0.6))
    seal(L)
    return L


# ---------------------------------------------------------------- C · Tobi-ishi
def mark_C():
    """The landscape printer: the stream runs left to right; to cross at the marked stone you
    look upstream. Five stepping stones: an upwind-biased stencil i-3 .. i+1."""
    L = Layer()
    amr(L, (32, 16, 64, 48), (40, 24, 56, 40))
    sd = 900
    # stream lines
    for k, (y0, amp, ph, w, x1) in enumerate(((11.0, 1.4, 0.0, 1.5, 61), (19.5, 1.8, 1.2, 1.9, 61), (44.0, 1.6, 2.1, 1.9, 50), (51.5, 1.2, 3.0, 1.4, 49))):
        pts = [(x, y0 + amp * math.sin(x / 7.5 + ph)) for x in np.linspace(3, x1, 30)]
        L.add("ink", stroke(pts, w, seed=sd + k, step=0.35, rise=0.06, tail=0.3, tail_min=0.12, bristles=5, kasure=0.35))
    # short current lines threading between stones
    xs = [11.0, 21.0, 31.0, 41.5, 52.0]
    yS = 31.5
    for k, (y0, xa, xb) in enumerate(((26.5, 4, 58), (36.5, 4, 56))):
        pts = [(x, y0 + 0.9 * math.sin(x / 5 + k)) for x in np.linspace(xa, xb, 24)]
        L.add("ink", stroke(pts, 0.6, seed=sd + 20 + k, step=0.35, rise=0.05, tail=0.3, tail_min=0.2, bristles=3, kasure=0.5))
    # eddies downstream of each stone
    sizes = [(3.3, 2.3), (3.6, 2.5), (3.4, 2.3), (4.0, 2.8), (3.0, 2.1)]
    for i, (x, (rx, ry)) in enumerate(zip(xs, sizes)):
        for j, rr in enumerate((rx + 1.6, rx + 3.0)):
            pts = arc(x + 0.8, yS, rr, -0.85, 0.85, 16)
            pts = [(px, yS + (py - yS) * 0.8) for px, py in pts]
            L.add("acc" if i == 3 else "ink", stroke(pts, 0.55, seed=sd + 40 + i * 3 + j, step=0.25, tail_min=0.2, wobble=0.1))
    for i, (x, (rx, ry)) in enumerate(zip(xs, sizes)):
        role = "acc" if i == 3 else "ink"
        L.add(role, dot(x, yS, rx, seed=sd + 60 + i, squash=ry / rx, angle=random.Random(i).uniform(-0.2, 0.2), rough=0.16))
    seal(L)
    return L


# ---------------------------------------------------------------- D · Ensō and tangent
def mark_D():
    """The Zen monk: the circle in one breath, the tangent touching once. Around the touch,
    the five neighbours that measure it."""
    L = Layer()
    amr(L, (0, 0, 32, 32), (8, 8, 24, 24))
    cx, cy, r = 33.0, 33.5, 19.0
    a0 = math.radians(-75)
    pts = arc(cx, cy, r, a0, a0 + math.radians(328), 90)
    prof = lambda t: (0.45 + 0.55 * np.clip(t / 0.08, 0, 1)) * (1 + 0.22 * np.sin(np.pi * t * 1.4)) * (1 - 0.7 * np.clip((t - 0.85) / 0.15, 0, 1))
    L.add("ink", stroke(pts, 5.2, seed=1200, step=0.35, profile=prof, wobble=0.06, bristles=7, kasure=0.32))
    th = math.radians(-142)
    px, py = cx + (r + 2.7) * math.cos(th), cy + (r + 2.7) * math.sin(th)
    tx, ty = -math.sin(th), math.cos(th)
    Lh = 21
    L.add("acc", stroke([(px - tx * Lh, py - ty * Lh), (px + tx * Lh * 0.82, py + ty * Lh * 0.82)], 2.0, seed=1300, step=0.3,
                        rise=0.05, tail=0.35, tail_min=0.15, press=0.35, bristles=3, kasure=0.15))
    for k in (-2, -1, 1, 2):
        a = th + k * math.radians(17)
        rr = r + 6.6
        L.add("ink", dot(cx + rr * math.cos(a), cy + rr * math.sin(a), 0.95 if abs(k) == 1 else 0.7, seed=1400 + k, squash=0.85, rough=0.15, tail=0.35, angle=a))
    seal(L)
    return L


# ---------------------------------------------------------------- E · Shōji
def mark_E():
    """The screen maker: a pine casts its shadow across five panes; the middle pane holds
    only the average of the light that falls on it. The bars between panes are the faces."""
    L = Layer()
    amr(L, (16, 16, 48, 48), (24, 16, 40, 32))
    x0, x1, y0, y1 = 5.5, 58.5, 12.5, 48.5
    xs = [x0 + (x1 - x0) * k / 5 for k in range(6)]
    # the average: a flat wash in the middle pane
    L.add("acc", poly([(xs[2] + 1.0, y0 + 1.1), (xs[3] - 1.0, y0 + 1.0), (xs[3] - 1.05, y1 - 1.0), (xs[2] + 1.05, y1 - 1.1)], seed=5, rough=0.12), 0.32)
    # pine shadow: a bough and needle tufts, faint
    bough = [(1, 22), (10, 19), (20, 21.5), (30, 19.5), (40, 23), (50, 21), (63, 24)]
    L.add("mesh", stroke(bough, 2.0, seed=1500, step=0.35, rise=0.05, tail=0.25, tail_min=0.2, bristles=4, kasure=0.3), 0.5)
    rng = random.Random(15)
    for k, bx in enumerate((7, 15, 24, 35, 45, 55)):
        by = 21 + rng.uniform(-1.5, 1.5)
        for j in range(5):
            a = math.radians(-150 + j * 30 + rng.uniform(-8, 8))
            ln = rng.uniform(3.5, 5.5)
            L.add("mesh", stroke([(bx, by), (bx + ln * math.cos(a), by + ln * math.sin(a) * 1.1 + 3.2)], 0.55, seed=1600 + k * 10 + j, step=0.2, tail_min=0.1), 0.5)
    # frame
    for (p, q, w) in (((x0 - 1.5, y0), (x1 + 1.5, y0), 2.4), ((x0 - 1.5, y1), (x1 - 3.0, y1), 2.4)):
        L.add("ink", stroke([p, q], w, seed=int(p[1] * 10), step=0.35, rise=0.04, tail=0.12, tail_min=0.6, press=0.3, bristles=5, kasure=0.12))
    for k, x in enumerate(xs):
        w = 2.1 if k in (0, 5) else 1.25
        ye = y1 + 1.0 if k < 5 else y1 - 3.5
        L.add("ink", stroke([(x, y0 - 0.8), (x + 0.15, ye)], w, seed=1700 + k, step=0.35, rise=0.05, tail=0.12, tail_min=0.6, press=0.2, bristles=3, kasure=0.08))
    # kumiko rail across the panes
    ym = (y0 + y1) / 2 + 2
    L.add("ink", stroke([(x0, ym), (x1, ym - 0.2)], 0.75, seed=1800, step=0.35, tail_min=0.7, wobble=0.1))
    seal(L)
    return L


# ---------------------------------------------------------------- F · Soroban
def mark_F():
    """The abacus master: nothing typed by hand, every weight counted bead by bead.
    The rods read 1 8 0 8 1, the weights of the 4th-order first derivative over twelve."""
    L = Layer()
    amr(L, (16, 16, 48, 48), (24, 16, 40, 32))
    x0, x1, y0, y1 = 6.0, 58.0, 11.0, 50.0
    beam = 22.5
    xs = [x0 + (x1 - x0) * (k + 0.5) / 5 for k in range(5)]
    digits = [1, 8, 0, 8, 1]
    # rods
    for k, x in enumerate(xs):
        L.add("acc" if k == 2 else "ink", stroke([(x, y0 + 0.5), (x + 0.1, y1 - 0.5)], 0.55 if k != 2 else 0.8, seed=1900 + k, step=0.35, tail_min=0.7, wobble=0.1))
    # frame
    L.add("ink", stroke([(x0 - 1.5, y0), (x1 + 1.2, y0 - 0.3)], 2.6, seed=1950, step=0.35, rise=0.04, tail=0.1, tail_min=0.6, press=0.3, bristles=5, kasure=0.12))
    L.add("ink", stroke([(x0 - 1.0, y1), (x1 - 4.5, y1 + 0.2)], 2.6, seed=1951, step=0.35, rise=0.04, tail=0.15, tail_min=0.4, press=0.3, bristles=5, kasure=0.25))
    L.add("ink", stroke([(x0, y0 - 1), (x0 + 0.2, y1 + 1)], 2.2, seed=1952, step=0.35, rise=0.04, tail=0.1, tail_min=0.6, bristles=4, kasure=0.1))
    L.add("ink", stroke([(x1, y0 - 1.2), (x1 - 0.2, y1 - 3.5)], 2.2, seed=1953, step=0.35, rise=0.04, tail=0.1, tail_min=0.5, bristles=4, kasure=0.15))
    L.add("ink", stroke([(x0 - 0.5, beam), (x1 + 0.5, beam - 0.2)], 1.5, seed=1954, step=0.35, rise=0.04, tail=0.1, tail_min=0.7, bristles=4, kasure=0.1))
    # unit dots on the beam, the centre one in colour
    for k, x in enumerate(xs):
        if k == 2:
            L.add("acc", dot(x, beam, 1.25, seed=1960, rough=0.1))
    bh, bw = 2.05, 3.7
    def bead(x, y, seed):
        return poly([(x - bw, y), (x - bw * 0.35, y - bh), (x + bw * 0.35, y - bh), (x + bw, y), (x + bw * 0.35, y + bh), (x - bw * 0.35, y + bh)], seed=seed, rough=0.15)
    for k, (x, dgt) in enumerate(zip(xs, digits)):
        role = "acc" if k == 2 else "ink"
        five = dgt >= 5
        ones = dgt - 5 * five
        yh = beam - 1.2 - bh if five else y0 + 1.6 + bh
        L.add(role, bead(x, yh, 2000 + k), 0.5 if role == "acc" else None)
        step = 2 * bh + 0.35
        for j in range(4):
            if j < ones:
                yb = beam + 1.1 + bh + j * step
            else:
                yb = y1 - 1.6 - bh - (3 - j) * step
            L.add(role, bead(x, yb, 2010 + k * 5 + j), 0.5 if role == "acc" else None)
    seal(L)
    return L


ARTISTS = [
    dict(id="A", file="Main.dc.html", fn=mark_A, role="The Garden Monk", name="Karesansui", kanji="枯山水",
         quote="Five places on the raked gravel. Four hold a stone; the fifth, the one being measured, holds only its ripples.",
         reads="The stones are the weights of the fourth-order first derivative, (+1, −8, 0, +8, −1)/12: big and small by magnitude, ink or pale by sign. The centre weight is zero, so the measured place is empty, and that emptiness is the colour. The rake lines bend around every stone: space is shaped by what sits in it.",
         prim="stone", palette="tokusa"),
    dict(id="B", file="Take.dc.html", fn=mark_B, role="The Bamboo Painter", name="Take no fushi", kanji="竹の節",
         quote="Bamboo grows in equal intervals. Each node knows only its neighbours, and from them the whole stalk can be measured.",
         reads="The stalk is the grid: five equal intervals of h, each node a cell centre. The pale stalk behind has its nodes halfway, at the faces. The sprig grows from the node being measured and leans the way the slope leans: one stroke in colour.",
         prim="node", palette="yamabuki"),
    dict(id="C", file="Tobiishi.dc.html", fn=mark_C, role="The Landscape Printer", name="Tobi-ishi", kanji="飛石",
         quote="The stream does not wait. To cross it at a stone, you look at the stones it comes from.",
         reads="The stream runs left to right; the stepping stones are an upwind-biased stencil, three upstream of the marked stone and one past it. Eddies trail each stone downstream, and the marked stone's eddies are in colour. Here space is the water you do not step on.",
         prim="flat", palette="ruri"),
    dict(id="D", file="Enso.dc.html", fn=mark_D, role="The Zen Monk", name="Ensō to sessen", kanji="円相と接線",
         quote="The circle is drawn in one breath. The tangent touches it only once, and in that touch is all its slope.",
         reads="The ensō is any smooth function; the coloured stroke is its derivative at one point, the only thing a derivative ever says. The four small dots outside are the neighbours that measure it. The empty middle of the circle is left alone, as the monks leave it.",
         prim="ring", palette="benibana"),
    dict(id="E", file="Shoji.dc.html", fn=mark_E, role="The Screen Maker", name="Shōji", kanji="障子",
         quote="A screen is mostly paper. What it shows is the light averaged over each pane, and the edges where one pane meets the next.",
         reads="Finite volumes, drawn as a sliding screen: five panes, a pine casting its shadow across them (the function), and the middle pane washed flat in colour (the cell average). The bars between panes are the faces where GIDEON reconstructs values.",
         prim="pane", palette="kikyo"),
    dict(id="F", file="Soroban.dc.html", fn=mark_F, role="The Abacus Master", name="Soroban", kanji="算盤",
         quote="Nothing typed by hand. Every weight is counted, bead by bead, and the count is exact.",
         reads="The rods read 1 · 8 · 0 · 8 · 1, the magnitudes of the fourth-order first-derivative weights over twelve, in real soroban notation. The middle rod counts zero and is in colour: the point being measured, with only the beam's unit dot on it. Exact rational arithmetic, counted rather than typed.",
         prim="bead", palette="hisui"),
]


def amr_dashed(L, r1, r2, seed=11):
    """Sumi-dashed AMR ground: every cell edge is its own brush dash, so the nodes are left
    as small gaps of bare paper. Three levels, three rhythms, three intensities:
      L0  h=16  two long wet dashes per edge, faint      (0.186)
      L1  h=8   one dash per edge, half the weight       (0.384)
      L2  h=4   short crisp ticks, the darkest           (0.62)"""
    rng = random.Random(seed)
    def dash(p, q, w, op, dry):
        prof = lambda t: (0.55 + 0.45 * np.clip(t / 0.18, 0, 1)) * (1 + 0.25 * np.exp(-((t - 0.1) / 0.08) ** 2)) * (1 - 0.75 * np.clip((t - 0.7) / 0.3, 0, 1))
        # a slight belly to each dash
        mx, my = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
        nx, ny = -(q[1] - p[1]), q[0] - p[0]
        ln = math.hypot(nx, ny) + 1e-9
        b = rng.uniform(-0.18, 0.18)
        m = (mx + nx / ln * b, my + ny / ln * b)
        L.add("mesh", stroke([p, m, q], w, seed=rng.randint(0, 10**6), step=0.3, profile=prof, wobble=0.12,
                             bristles=3 if dry else 1, kasure=dry), op)
    def edges(h, a0, b0, a1, b1, skip, per, gap, w, op, dry, horizontal):
        lines = [v for v in range((b0 if horizontal else a0), (b1 if horizontal else a1) + 1, h) if v not in skip]
        lo, hi = (a0, a1) if horizontal else (b0, b1)
        for v in lines:
            if v <= 0 or v >= 64:
                continue
            for s in range(lo, hi, h):
                seg = h / per
                for k in range(per):
                    u0 = s + k * seg + gap * rng.uniform(0.7, 1.3)
                    u1 = s + (k + 1) * seg - gap * rng.uniform(0.7, 1.3)
                    if rng.random() < 0.5:  # brush direction alternates like a hand would
                        u0, u1 = u1, u0
                    j = rng.uniform(-0.12, 0.12)
                    if horizontal:
                        dash((u0, v + j), (u1, v - j), w, op, dry)
                    else:
                        dash((v + j, u0), (v - j, u1), w, op, dry)
    for hz in (True, False):
        edges(16, 0, 0, 64, 64, set(), 2, 1.0, 0.7, AMR_OPS[0], 0.25, hz)
    a0, b0, a1, b1 = r1
    for hz in (True, False):
        skip = {v for v in range(0, 65, 16)}
        edges(8, a0, b0, a1, b1, skip, 1, 0.8, 0.6, AMR_OPS[1], 0.0, hz)
    a0, b0, a1, b1 = r2
    for hz in (True, False):
        skip = {v for v in range(0, 65, 8)}
        edges(4, a0, b0, a1, b1, skip, 1, 0.6, 0.5, AMR_OPS[2], 0.0, hz)


def mark_D_final():
    L = mark_D()
    L.parts = [q for q in L.parts if not (q[0] == "mesh" and q[1] in AMR_OPS)]
    G = Layer()
    amr_dashed(G, (0, 0, 32, 32), (8, 8, 24, 24))
    L.parts = G.parts + L.parts
    return L
