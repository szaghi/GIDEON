"""Nine-icon pack, 24x24, drawn with each artist's primitive."""
import math
import numpy as np
from brush import stroke, dot, poly, arc, Layer


def P(L, kind, x, y, s, role, seed, op=None):
    if kind == "stone":
        L.add(role, dot(x, y, s, seed=seed, squash=0.78, angle=0.3, rough=0.12), op)
    elif kind == "node":
        L.add(role, stroke([(x - s * 1.25, y + s * 0.2), (x, y - s * 0.15), (x + s * 1.25, y + s * 0.2)], s * 0.95, seed=seed, step=0.15, rise=0.15, tail=0.3, tail_min=0.3, press=0.3), op)
    elif kind == "flat":
        L.add(role, dot(x, y, s * 1.2, seed=seed, squash=0.62, angle=-0.1, rough=0.17), op)
    elif kind == "ring":
        a0 = (seed % 7) * 0.9
        L.add(role, stroke(arc(x, y, s * 0.82, a0, a0 + 2 * math.pi * 0.86, 24), s * 0.55, seed=seed, step=0.12, rise=0.1, tail=0.3, tail_min=0.25), op)
    elif kind == "pane":
        h = s * 0.95
        L.add(role, poly([(x - h, y - h), (x + h, y - h * 1.05), (x + h * 1.03, y + h), (x - h * 0.98, y + h)], seed=seed, rough=0.12), op)
    elif kind == "bead":
        bw, bh = s * 1.35, s * 0.8
        L.add(role, poly([(x - bw, y), (x - bw * 0.35, y - bh), (x + bw * 0.35, y - bh), (x + bw, y), (x + bw * 0.35, y + bh), (x - bw * 0.35, y + bh)], seed=seed, rough=0.08), op)


def S(L, pts, w, role, seed, op=None, **kw):
    kw.setdefault("step", 0.18)
    kw.setdefault("tail_min", 0.3)
    L.add(role, stroke(pts, w, seed=seed, **kw), op)


ICONS = ["Exact coefficients", "Derivatives d1–d6", "Face values", "Cell averages",
         "Centred and upwind", "Device-callable kernels", "Grad, div, curl, Laplacian", "Verified orders", "Multi-licensed"]


def icon(kind, i, wm=1.0):
    L = Layer()
    sd = 31 * i + 7
    if i == 0:  # exact: p / q
        S(L, [(4.5, 12.2), (19.5, 11.8)], 2.0 * wm, "ink", sd, rise=0.05, press=0.3)
        P(L, kind, 12, 6.0, 2.0, "acc", sd + 1)
        P(L, kind, 12, 18.0, 2.0, "ink", sd + 2)
    elif i == 1:  # derivative: a hill and its tangent
        f = lambda x: 9.5 + 0.085 * (x - 13) ** 2
        S(L, [(x, f(x)) for x in np.linspace(3, 22, 12)], 1.7 * wm, "ink", sd)
        x0 = 7.5; y0 = f(x0); k = 0.17 * (x0 - 13)
        S(L, [(1.5, y0 + k * (1.5 - x0) - 1.1), (13.5, y0 + k * (13.5 - x0) - 1.1)], 1.3 * wm, "acc", sd + 1, rise=0.05, tail_min=0.2)
        P(L, kind, x0 + 0.3, y0 + 1.0, 1.25, "ink", sd + 2)
    elif i == 2:  # face between two cells
        P(L, kind, 6.0, 12, 2.1, "ink", sd)
        P(L, kind, 18.0, 12, 2.1, "ink", sd + 1)
        S(L, [(12.1, 3.5), (11.9, 20.5)], 1.9 * wm, "acc", sd + 2, rise=0.05, press=0.3)
    elif i == 3:  # cell average: brackets, a wave, its mean
        S(L, [(6.5, 4.5), (4.2, 12), (6.5, 19.5)], 1.6 * wm, "ink", sd, tail_min=0.4)
        S(L, [(17.5, 4.5), (19.8, 12), (17.5, 19.5)], 1.6 * wm, "ink", sd + 1, tail_min=0.4)
        S(L, [(8, 9.5), (10, 7.8), (12, 9.8), (14, 11.0), (16, 9.2)], 0.9 * wm, "mesh", sd + 2, op=0.5 if False else None)
        S(L, [(8, 15.0), (16, 15.0)], 1.7 * wm, "acc", sd + 3, rise=0.05, tail_min=0.5)
    elif i == 4:  # centred and upwind rows
        for k, x in enumerate((3.5, 8, 12.5 - 0.5, 16, 20.5)):
            P(L, kind, [3.5, 7.75, 12, 16.25, 20.5][k], 7.5, 1.25, "acc" if k == 2 else "ink", sd + k)
        for k, x in enumerate((3.5, 7.75, 12, 16.25)):
            P(L, kind, x, 16.5, 1.25, "acc" if k == 3 else "ink", sd + 10 + k)
        S(L, [(3, 20.8), (14.5, 20.8)], 0.9 * wm, "mesh", sd + 20, tail_min=0.15)
    elif i == 5:  # device kernel: a chip with pins
        for a, b in (((6, 6), (18, 6)), ((18, 6), (18, 18)), ((18, 18), (6, 18)), ((6, 18), (6, 6))):
            S(L, [a, b], 1.5 * wm, "ink", sd + int(a[0] + b[1]), rise=0.05, tail_min=0.6)
        for y in (9, 12, 15):
            S(L, [(2, y), (4.6, y)], 1.0 * wm, "ink", sd + y, tail_min=0.5)
            S(L, [(19.4, y), (22, y)], 1.0 * wm, "ink", sd + y + 50, tail_min=0.5)
        P(L, kind, 12, 12, 2.0, "acc", sd + 99)
    elif i == 6:  # nabla
        S(L, [(3.5, 5.0), (20.5, 5.0)], 1.8 * wm, "ink", sd, rise=0.05, press=0.3, tail_min=0.5)
        S(L, [(20.5, 5.2), (12.2, 20.5)], 1.6 * wm, "ink", sd + 1, tail_min=0.4)
        S(L, [(3.6, 5.2), (11.8, 20.5)], 1.3 * wm, "ink", sd + 2, tail_min=0.4)
        P(L, kind, 12, 10.2, 1.6, "acc", sd + 3)
    elif i == 7:  # verified orders: staircase and a tick
        stair = [(2.5, 8.5), (2.5, 13), (7, 13), (7, 17), (11.5, 17), (11.5, 21), (17.5, 21)]
        for q in range(len(stair) - 1):
            S(L, [stair[q], stair[q + 1]], 1.2 * wm, "ink", sd + 40 + q, rise=0.05, tail=0.15, tail_min=0.55)
        S(L, [(10.5, 8.5), (14, 12), (21.5, 3.0)], 2.1 * wm, "acc", sd + 1, rise=0.1, tail=0.4, tail_min=0.2, press=0.3)
    elif i == 8:  # multi-licensed: seal with four cut primitives
        L.add("acc", poly([(4.2, 4.3), (19.8, 4.0), (20.0, 19.9), (4.0, 20.1)], seed=sd, rough=0.15))
        for k, (x, y) in enumerate(((8.6, 8.6), (15.4, 8.6), (8.6, 15.4), (15.4, 15.4))):
            P(L, kind, x, y, 1.9, "ground", sd + k)
    return L
