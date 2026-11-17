import os, json, math, shutil, subprocess
from brush import stroke, dot, arc, Layer, d_of
from marks import mark_D_final, seal, AMR_OPS
from icons import icon, ICONS
from identity import W, N, CSS, cr

OUT = os.path.join(os.environ.get("GIDEON_BRAND_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")), "gideon-identity")  # was /home/claude/bundle/gideon-identity
shutil.rmtree(OUT, ignore_errors=True)
PUB = f"{OUT}/docs/public"; THEME = f"{OUT}/docs/.vitepress/theme"
os.makedirs(f"{PUB}/icons", exist_ok=True); os.makedirs(THEME, exist_ok=True)

def cols(p): return dict(ink=p["ink"], acc=p["acc"], mesh=p["mesh"], ground=p["ground"])

def svgfile(L, p, vb, label, rect=None, prec=2):
    s = L.svg(cols(p), viewBox=vb, size=0, label=label, prec=prec)
    s = s.replace(' width="0" height="0"', '')
    if rect:
        x, y, w, r = rect
        s = s.replace(f'aria-label="{label}">', f'aria-label="{label}"><rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{w:.2f}" rx="{r:.2f}" fill="{p["raised"]}"/>', 1)
    return s

def small_mark():
    """Navbar and favicon: ensō and tangent only, strokes thickened for 16-32 px."""
    L = Layer()
    cx, cy, r = 33.0, 33.5, 19.0
    a0 = math.radians(-75)
    prof = lambda t: (0.6 + 0.4 * min(1, t / 0.08) if False else 1) * 1.0
    import numpy as np
    prof = lambda t: (0.55 + 0.45 * np.clip(t / 0.08, 0, 1)) * (1 + 0.15 * np.sin(np.pi * t * 1.4)) * (1 - 0.55 * np.clip((t - 0.88) / 0.12, 0, 1))
    L.add("ink", stroke(arc(cx, cy, r, a0, a0 + math.radians(328), 90), 6.6, seed=1200, step=0.35, profile=prof, wobble=0.04))
    th = math.radians(-142)
    px, py = cx + (r + 3.6) * math.cos(th), cy + (r + 3.6) * math.sin(th)
    tx, ty = -math.sin(th), math.cos(th)
    Lh = 20
    L.add("acc", stroke([(px - tx * Lh, py - ty * Lh), (px + tx * Lh * 0.8, py + ty * Lh * 0.8)], 3.6, seed=1300, step=0.3, rise=0.05, tail=0.35, tail_min=0.3, press=0.3))
    return L

M = mark_D_final()
S = small_mark()
x0, y0, x1, y1 = S.bbox(("ink", "acc")); s = max(x1 - x0, y1 - y0) + 1.0
cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
NAV_VB = f"{cx - s/2:.2f} {cy - s/2:.2f} {s:.2f} {s:.2f}"
T = s * 1.18
FAV_VB = f"{cx - T/2:.2f} {cy - T/2:.2f} {T:.2f} {T:.2f}"
files = {}
for suf, p in (("", W), ("-dark", N)):
    files[f"logo{suf}.svg"] = svgfile(M, p, "0 0 64 64", "GIDEON")
    files[f"logo-nav{suf}.svg"] = svgfile(S, p, NAV_VB, "GIDEON")
    files[f"favicon{suf}.svg"] = svgfile(S, p, FAV_VB, "GIDEON", rect=(cx - T/2, cy - T/2, T, T * 0.19))
NAMES = ["exact", "derivatives", "faces", "averages", "stencils", "kernels", "nabla", "verified", "licence"]
for i, nm in enumerate(NAMES):
    L = icon("ring", i)
    for suf, p in (("", W), ("-dark", N)):
        files[f"icons/{nm}{suf}.svg"] = svgfile(L, p, "0 0 24 24", ICONS[i], prec=2)
for fn, t in files.items():
    open(f"{PUB}/{fn}", "w").write(t)
for px in (16, 32):
    subprocess.run(["rsvg-convert", "-w", str(px), "-h", str(px), f"{PUB}/favicon.svg", "-o", f"{PUB}/favicon-{px}.png"], check=True)
subprocess.run(["rsvg-convert", "-w", "180", "-h", "180", "-b", W["raised"], f"{PUB}/favicon.svg", "-o", f"{PUB}/apple-touch-icon.png"], check=True)

open(f"{THEME}/custom.css", "w").write(CSS)
open(f"{THEME}/index.ts", "w").write("// GIDEON docs theme: VitePress default theme + the Sumi family tokens in the Kogane ink.\nimport DefaultTheme from 'vitepress/theme'\nimport './custom.css'\n\nexport default DefaultTheme\n")

keys = ["ground", "alt", "raised", "line", "ink", "muted", "mesh", "acc", "link", "pressed", "on", "focus"]
tok = {"name": "GIDEON · Ensō to sessen · Kogane 黄金", "family": "Sumi family",
       "fonts": {"display": "Shippori Mincho 800", "sans": "Inter", "mono": "JetBrains Mono"},
       "amr_ground": {"style": "sumi-dashed, nodes left bare", "levels": [{"h": 16, "opacity": AMR_OPS[0]}, {"h": 8, "opacity": AMR_OPS[1]}, {"h": 4, "opacity": AMR_OPS[2]}]},
       "light": {k: W[k] for k in keys}, "dark": {k: N[k] for k in keys},
       "contrast": {m: {k: round(cr(p[k] if k != "on" else p["on"], p["acc"] if k == "on" else p["ground"]), 2) for k in ("ink", "muted", "mesh", "acc", "link", "on", "focus")} for m, p in (("light", W), ("dark", N))}}
open(f"{OUT}/tokens.json", "w").write(json.dumps(tok, indent=2, ensure_ascii=False))
os.makedirs(f"{OUT}/generator", exist_ok=True)
for f in ("brush.py", "marks.py", "icons.py", "pal.py", "final.py", "identity.py", "bundle.py"):
    shutil.copy(f, f"{OUT}/generator/{f}")
print(sorted(os.listdir(PUB)), NAV_VB, FAV_VB)
