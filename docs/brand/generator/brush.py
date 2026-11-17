"""Sumi brush primitives: every mark is a set of filled polygons."""
import math, random
import numpy as np


def _smooth(a, b, x):
    t = np.clip((x - a) / max(b - a, 1e-9), 0, 1)
    return t * t * (3 - 2 * t)


def catmull(pts, step):
    """Resample a polyline through pts with a Catmull-Rom spline, ~step apart."""
    P = np.array(pts, float)
    if len(P) == 2:
        L = np.linalg.norm(P[1] - P[0])
        n = max(int(L / step), 2)
        t = np.linspace(0, 1, n)[:, None]
        return P[0] * (1 - t) + P[1] * t
    P = np.vstack([2 * P[0] - P[1], P, 2 * P[-1] - P[-2]])
    out = []
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        L = np.linalg.norm(p2 - p1)
        n = max(int(L / step), 1)
        for k in range(n):
            t = k / n
            t2, t3 = t * t, t * t * t
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 + (-p0 + 3 * p1 - 3 * p2 + p3) * t3))
    out.append(P[-2])
    return np.array(out)


def _normals(C):
    d = np.gradient(C, axis=0)
    d /= np.linalg.norm(d, axis=1)[:, None] + 1e-12
    return np.stack([-d[:, 1], d[:, 0]], 1)


def _noise(n, rng, amp, waves=3):
    t = np.linspace(0, 1, n)
    s = np.zeros(n)
    for k in range(waves):
        f = rng.uniform(1.5, 5.5) * (k + 1)
        s += math.sin(1) * 0 + np.sin(2 * math.pi * f * t + rng.uniform(0, 6.3)) / (k + 1)
    return 1 + amp * s / 1.6


def stroke(pts, w, seed=1, step=0.35, rise=0.12, tail=0.25, tail_min=0.18, press=0.15,
           wobble=0.08, bristles=1, kasure=0.0, profile=None, closed=False):
    """A brushed stroke along pts. Returns a list of polygons (arrays Nx2)."""
    rng = random.Random(seed)
    C = catmull(pts, step)
    n = len(C)
    if n < 3:
        return []
    N = _normals(C)
    t = np.linspace(0, 1, n)
    if profile is not None:
        W = w * profile(t)
    else:
        W = w * (0.4 + 0.6 * _smooth(0, rise, t)) * (1 - (1 - tail_min) * _smooth(1 - tail, 1, t))
        W *= 1 + press * np.exp(-((t - 0.07) / 0.06) ** 2)
    W *= _noise(n, rng, wobble)
    # slight lateral drift of the centreline
    drift = (_noise(n, rng, 1.0, 2) - 1) * w * 0.06
    C = C + N * drift[:, None]
    polys = []
    if bristles <= 1:
        L = C + N * (W / 2)[:, None]
        R = C - N * (W / 2)[:, None]
        polys.append(np.vstack([L, R[::-1]]))
        return polys
    k = bristles
    # solid core so the stroke never reads as stripes
    core = 0.62
    Lc = C + N * (W * core / 2)[:, None]
    Rc = C - N * (W * core / 2)[:, None]
    tcut = 1 - kasure * 0.55 if kasure > 0 else 1
    mc = max(int(tcut * (n - 1)), 2)
    polys.append(np.vstack([Lc[:mc], Rc[:mc][::-1]]))
    for j in range(k):
        u0 = -1 + 2 * j / k - 0.04
        u1 = -1 + 2 * (j + 1) / k + 0.04
        if kasure > 0 and rng.random() < 0.7:
            cut = 1 - kasure * (rng.random() ** 0.6)
        else:
            cut = 1 - rng.uniform(0, 0.04)
        m = max(int(cut * (n - 1)), 2)
        tt = np.linspace(0, 1, m)
        shrink = 1 - 0.85 * _smooth(0.9, 1, tt)
        mid = (u0 + u1) / 2
        half = (u1 - u0) / 2 * shrink
        a = (mid - half) * (W[:m] / 2)
        b = (mid + half) * (W[:m] / 2)
        polys.append(np.vstack([C[:m] + N[:m] * b[:, None], (C[:m] + N[:m] * a[:, None])[::-1]]))
        if kasure > 0 and cut < 0.95 and rng.random() < 0.45:
            s0 = int((cut + 0.025) * (n - 1))
            s1 = min(int((cut + 0.025 + rng.uniform(0.03, 0.1)) * (n - 1)), n - 1)
            if s1 - s0 > 3:
                sl = slice(s0, s1)
                hh = (u1 - u0) / 2 * 0.45
                polys.append(np.vstack([C[sl] + N[sl] * ((mid + hh) * W[sl] / 2)[:, None], (C[sl] + N[sl] * ((mid - hh) * W[sl] / 2)[:, None])[::-1]]))
    return polys


def dot(x, y, r, seed=1, squash=0.85, angle=0.0, rough=0.08, tail=0.0, n=28):
    """A ten (brush dot) or stone: a rough blob."""
    rng = random.Random(seed)
    ph = [rng.uniform(0, 6.3) for _ in range(3)]
    pts = []
    ca, sa = math.cos(angle), math.sin(angle)
    for i in range(n):
        th = 2 * math.pi * i / n
        rr = r * (1 + rough * (math.sin(2 * th + ph[0]) * 0.6 + math.sin(3 * th + ph[1]) * 0.4 + math.sin(5 * th + ph[2]) * 0.25))
        # tail: elongate towards angle 0
        rr *= 1 + tail * max(0, math.cos(th)) ** 3
        px, py = rr * math.cos(th), rr * squash * math.sin(th)
        pts.append((x + px * ca - py * sa, y + px * sa + py * ca))
    return [np.array(pts)]


def poly(points, seed=1, rough=0.0):
    rng = random.Random(seed)
    return [np.array([(x + rng.uniform(-rough, rough), y + rng.uniform(-rough, rough)) for x, y in points])]


def arc(cx, cy, r, a0, a1, n=40):
    return [(cx + r * math.cos(a0 + (a1 - a0) * i / (n - 1)), cy + r * math.sin(a0 + (a1 - a0) * i / (n - 1))) for i in range(n)]


def rdp(P, eps):
    P = np.asarray(P)
    if len(P) < 4:
        return P
    keep = np.zeros(len(P), bool); keep[0] = keep[-1] = True
    stack = [(0, len(P) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        a, b = P[i], P[j]
        ab = b - a; L = np.hypot(*ab)
        seg = P[i + 1:j] - a
        if L < 1e-9:
            d = np.hypot(seg[:, 0], seg[:, 1])
        else:
            d = np.abs(ab[0] * seg[:, 1] - ab[1] * seg[:, 0]) / L
        k = int(np.argmax(d))
        if d[k] > eps:
            keep[i + 1 + k] = True
            stack += [(i, i + 1 + k), (i + 1 + k, j)]
    return P[keep]


EPS = [0.035]

def d_of(polys, prec=2):
    f = "{:.%df}" % prec
    out = []
    for P in polys:
        P = rdp(P, EPS[0])
        if len(P) < 3:
            continue
        s = "M" + f.format(P[0][0]) + " " + f.format(P[0][1]) + "L" + " ".join(f.format(x) + " " + f.format(y) for x, y in P[1:]) + "Z"
        out.append(s)
    return " ".join(out)


class Layer:
    """Collects polygons by colour role: 'ink', 'acc', 'mesh', 'ground', and AMR levels."""

    def __init__(self):
        self.parts = []  # (role, opacity, polys, evenodd)

    def add(self, role, polys, op=None, evenodd=False):
        self.parts.append((role, op, polys, evenodd))

    def svg(self, colors, viewBox="0 0 64 64", size=200, label="", prec=2, roles=None):
        groups = {}
        order = []
        for role, op, polys, eo in self.parts:
            if roles is not None and role not in roles:
                continue
            key = (role, op, eo)
            if key not in groups:
                groups[key] = []
                order.append(key)
            groups[key].extend(polys)
        body = []
        for key in order:
            role, op, eo = key
            col = colors[role]
            attrs = f'fill="{col}"'
            if op is not None:
                attrs += f' opacity="{op}"'
            if eo:
                attrs += ' fill-rule="evenodd"'
            body.append(f'<path {attrs} d="{d_of(groups[key], prec)}"/>')
        lab = f'role="img" aria-label="{label}"' if label else 'aria-hidden="true"'
        return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{viewBox}" width="{size}" height="{size}" {lab}>' + "".join(body) + "</svg>"

    def bbox(self, roles):
        xs, ys = [], []
        for role, op, polys, eo in self.parts:
            if role in roles:
                for P in polys:
                    xs.extend(P[:, 0]); ys.extend(P[:, 1])
        return min(xs), min(ys), max(xs), max(ys)
