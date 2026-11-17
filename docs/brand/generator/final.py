import json, os, re
from brush import Layer
from marks import mark_D, mark_D_final, AMR_OPS
from icons import icon, ICONS
from pal import cr, fit, shift, WASHI, NIGHT

ROOT = os.path.join(os.environ.get("GIDEON_BRAND_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")), "boards")  # was /home/claude/gideon/project
FONTS = '<link rel="preconnect" href="https://fonts.googleapis.com"/><link href="https://fonts.googleapis.com/css2?family=Shippori+Mincho:wght@800&amp;family=Inter:wght@400;500;600&amp;family=JetBrains+Mono&amp;display=swap" rel="stylesheet"/>'
DISP = "font-family: 'Shippori Mincho', 'Hiragino Mincho ProN', serif"
MONO = "font-family: 'JetBrains Mono', ui-monospace, monospace"
TAG = "Generated Interpolations and Derivatives, Exact to Order N"

RAW = [
 ("benibana", "Benibana", "紅花", "Safflower crimson, the red the ensō was first shown in.", {}, {}, "#b0243f", "#ff7f96"),
 ("urushi", "Urushi", "漆", "Vermilion lacquer. At night the ground itself turns to black lacquer, warm and deep.", {},
   dict(ground="#160e0b", alt="#110a08", raised="#21160f", line="#36251d", ink="#f2e6dc", muted="#c4ae9f", mesh="#9c8678"), "#bf3a1b", "#ff6b3d"),
 ("kogane", "Kogane", "黄金", "Gold. The tangent as a line of gold leaf; at night, gold on smoked bronze.", {},
   dict(ground="#11100c", alt="#0c0b08", raised="#1b1913", line="#2e2a1f", ink="#efe9da", muted="#bab29e", mesh="#948c79"), "#8a6300", "#f2c14e"),
 ("seiboku", "Seiboku", "青墨", "Blue-black ink. Almost monochrome: the tangent is only a cooler, glossier ink than the circle.", {}, {}, "#2b4560", "#b4cbe2"),
 ("aigami", "Aigami", "藍紙", "Indigo paper. Light mode on pale indigo-washed paper; at night, deep indigo with a persimmon tangent.",
   dict(ground="#e9eef4", alt="#dde4ee", raised="#f2f5f9", line="#cdd6e3", ink="#141c2b", muted="#4a5568", mesh="#66718a"),
   dict(ground="#0b1633", alt="#08112a", raised="#12204a", line="#24335f", ink="#eef2ff", muted="#aab6d6", mesh="#8090b5"), "#233f96", "#ff9a5c"),
 ("uguisu", "Uguisu", "鶯", "Warbler olive, the green-brown of the bush warbler in early spring.", {}, {}, "#56661a", "#c9da6c"),
 ("tsutsuji", "Neon Tsutsuji", "躑躅", "Electric azalea. A magenta that glows on a violet-black night.", {},
   dict(ground="#0f0b13", alt="#0b080e", raised="#18121e", line="#2b2134", ink="#f1eaf5", muted="#b8abc2", mesh="#9284a0"), "#ad1680", "#ff4fd2"),
 ("mizu", "Neon Mizu", "水", "Electric water. Cyan light, like a GPU running in a dark room.", {},
   dict(ground="#071114", alt="#050d10", raised="#0d1b20", line="#1b2e34", ink="#e6f6f8", muted="#9fbac0", mesh="#7b989e"), "#00707d", "#3ef0ff"),
]
PALS = []
for key, name, kj, desc, wo, no, a, b in RAW:
    w = dict(WASHI, **wo); n = dict(NIGHT, **no)
    w["acc"] = fit(a, w["ground"], 4.8, -1); w["link"] = fit(shift(w["acc"], -0.03), w["ground"], 5.5, -1)
    n["acc"] = fit(b, n["ground"], 7.0, +1); n["link"] = fit(shift(n["acc"], 0.04), n["ground"], 8.0, +1)
    for p in (w, n):
        p["cr"] = round(cr(p["acc"], p["ground"]), 1); p["crl"] = round(cr(p["link"], p["ground"]), 1)
        p["crink"] = round(cr(p["ink"], p["ground"]), 1); p["crmut"] = round(cr(p["muted"], p["ground"]), 1)
    PALS.append(dict(key=key, name=name, kanji=kj, desc=desc, washi=w, night=n))

M = mark_D_final()
OLD = mark_D()
F = Layer(); F.parts = [q for q in M.parts if q[1] not in AMR_OPS and not q[3]]
x0, y0, x1, y1 = F.bbox(("ink", "acc", "mesh")); s = max(x1 - x0, y1 - y0) + 1.6
FVB = f"{(x0+x1)/2 - s/2:.2f} {(y0+y1)/2 - s/2:.2f} {s:.2f} {s:.2f}"
ICN = [icon("ring", i) for i in range(9)]

def cols(p): return dict(ink=p["ink"], acc=p["acc"], mesh=p["mesh"], ground=p["ground"])
def mark(p, size, label, L=None): return (L or M).svg(cols(p), size=size, label=label, prec=1)
def favs(p, sizes=((64, 12, 54), (32, 6, 27), (16, 3, 14))):
    return '<div style="display: flex; gap: 8px; align-items: flex-end">' + "".join(
        f'<div style="width: {px}px; height: {px}px; border-radius: {r}px; background: {p["raised"]}; flex-shrink: 0; display: flex; align-items: center; justify-content: center; border: 1px solid {p["line"]}">{F.svg(cols(p), viewBox=FVB, size=inner, label=f"favicon {px}px", prec=1)}</div>'
        for px, r, inner in sizes) + "</div>"
def sw(p):
    return '<div style="display: flex; gap: 4px">' + "".join(
        f'<span style="display: inline-block; width: 18px; height: 14px; border-radius: 4px; background: {p[k]}; border: 1px solid {p["line"]}"></span>'
        for k in ("ground", "raised", "ink", "muted", "mesh", "acc", "link")) + "</div>"

def page(title, w, h, body, bg="#f3ede1", link="#1d1b18"):
    return f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<script src="./support.js"></script>
</head>
<body>
<x-dc>
<helmet>
{FONTS}
<style>
body{{margin:0;background:{bg}}}
a{{color:{link}}}a:hover{{color:#1d1b18}}
</style>
</helmet>
{body}
</x-dc>
<script type="text/x-dc" data-dc-script data-props='{{"$preview":{{"width":{w},"height":{h}}}}}'>
class Component extends DCLogic {{
renderVals() {{
return {{}};
}}
}}
</script>
</body>
</html>
'''

def ink_board(i, P):
    w, n = P["washi"], P["night"]
    cap = "margin: 0; font-size: 12px; line-height: 16px; font-weight: 500; letter-spacing: 0.02em"
    iw = "".join(f'<div style="display: flex; flex-direction: column; gap: 5px; align-items: center; width: 70px">{ICN[k].svg(cols(w), viewBox="0 0 24 24", size=32, label=ICONS[k], prec=1)}<span style="font-size: 10.5px; line-height: 13px; color: {w["muted"]}; text-align: center">{ICONS[k]}</span></div>' for k in range(9))
    inn = "".join(ICN[k].svg(cols(n), viewBox="0 0 24 24", size=26, label=ICONS[k], prec=1) for k in range(9))
    body = f'''<div style="width: 440px; height: 1000px; display: flex; flex-direction: column; font-family: Inter, ui-sans-serif, system-ui, sans-serif">
<div style="background: {w["ground"]}; color: {w["ink"]}; padding: 24px 28px 20px; display: flex; flex-direction: column; gap: 12px; flex-grow: 1">
<div style="display: flex; justify-content: space-between; align-items: baseline"><p style="{cap}; color: {w["muted"]}">Ink {i + 1} of 8 · Ensō to sessen</p><p style="{cap}; color: {w["link"]}">{P["name"]} {P["kanji"]}</p></div>
<h1 style="margin: 0; {DISP}; font-size: 28px; line-height: 32px; font-weight: 800; letter-spacing: 0.05em">{P["name"]} <span style="font-size: 18px; letter-spacing: 0.1em; color: {w["muted"]}">{P["kanji"]}</span></h1>
<p style="margin: 0; font-size: 13px; line-height: 19px; font-style: italic; color: {w["muted"]}">{P["desc"]}</p>
<div style="width: 384px; height: 384px">{mark(w, 384, "GIDEON logo, Ensō to sessen, " + P["name"] + ", washi")}</div>
<div style="display: flex; align-items: flex-end; justify-content: space-between"><div style="display: flex; flex-direction: column; gap: 6px"><span style="{DISP}; font-size: 36px; line-height: 36px; font-weight: 800; letter-spacing: 0.04em">GIDEON</span><span style="font-size: 11px; line-height: 15px; font-weight: 500; letter-spacing: 0.02em; color: {w["muted"]}">{TAG}</span></div></div>
<div style="display: flex; align-items: center; justify-content: space-between">{favs(w)}<div style="display: flex; flex-direction: column; gap: 6px; align-items: flex-end">{sw(w)}<span style="{MONO}; font-size: 10.5px; color: {w["muted"]}">{w["acc"]} {w["cr"]}:1 · link {w["crl"]}:1</span></div></div>
<div style="display: flex; flex-wrap: wrap; gap: 10px 4px">{iw}</div>
</div>
<div style="background: {n["ground"]}; color: {n["ink"]}; padding: 20px 28px 22px; display: flex; flex-direction: column; gap: 14px">
<div style="display: flex; gap: 16px; align-items: center">
<div style="width: 128px; height: 128px; flex-shrink: 0">{mark(n, 128, "GIDEON logo at night, " + P["name"])}</div>
<div style="display: flex; flex-direction: column; gap: 8px">
<span style="{DISP}; font-size: 30px; line-height: 30px; font-weight: 800; letter-spacing: 0.05em">GIDEON</span>
{favs(n)}
<span style="{MONO}; font-size: 10.5px; color: {n["muted"]}">{n["acc"]} {n["cr"]}:1 · link {n["crl"]}:1</span>
</div></div>
<div style="display: flex; justify-content: space-between">{inn}</div>
</div>
</div>'''
    return page(f"GIDEON · {P['name']}", 440, 1000, body, bg=w["ground"], link=w["link"])

def level_layer(op):
    L = Layer(); L.parts = [q for q in M.parts if q[0] == "mesh" and q[1] == op]; return L

def ground_board():
    w = PALS[0]["washi"]
    def card(svg, title, note):
        return f'<div style="display: flex; flex-direction: column; gap: 8px">{svg}<span style="font-size: 13px; font-weight: 600; color: #1d1b18">{title}</span><span style="font-size: 12px; line-height: 17px; color: #5f5a52">{note}</span></div>'
    lv = [(AMR_OPS[0], "L0 · h = 16 · 0.186", "Two long wet dashes per cell edge, a little dry at the tail. The whole square."),
          (AMR_OPS[1], "L1 · h = 8 · 0.384", "One dash per edge, firmer. Only where the tangent is."),
          (AMR_OPS[2], "L2 · h = 4 · 0.62", "Short crisp ticks, the darkest, around the point of contact.")]
    small = "".join(card(f'<div style="width: 200px; height: 200px; background: #f7f2e8; border: 1px solid #ddd4c3; border-radius: 8px">{level_layer(op).svg(cols(w), size=200, label=t, prec=1)}</div>', t, nt) for op, t, nt in lv)
    body = f'''<div style="width: 960px; height: 1000px; box-sizing: border-box; padding: 32px; background: #f3ede1; font-family: Inter, ui-sans-serif, system-ui, sans-serif; display: flex; flex-direction: column; gap: 22px">
<div style="display: flex; flex-direction: column; gap: 6px"><span style="{DISP}; font-size: 30px; line-height: 34px; font-weight: 800; letter-spacing: 0.04em; color: #1d1b18">The AMR ground, sumi-dashed</span>
<span style="font-size: 13px; line-height: 19px; color: #5f5a52; max-width: 860px">Every cell edge is now its own brush dash, and the grid nodes are left as small gaps of bare paper: the points of a stencil, with the space between them painted. Each level has its own rhythm and its own intensity, and the refinement still follows the tangent.</span></div>
<div style="display: flex; gap: 32px">
{card(f'<div style="width: 420px; height: 420px">{OLD.svg(cols(w), size=420, label="Round 1 ground", prec=1)}</div>', "Before · continuous hairlines", "Round 1: three levels of thin continuous lines.")}
{card(f'<div style="width: 420px; height: 420px">{mark(w, 420, "New sumi-dashed ground")}</div>', "After · sumi-dashed", "Three levels of dashes; nodes left open.")}
</div>
<div style="display: flex; gap: 32px">{small}</div>
</div>'''
    return page("GIDEON · AMR ground", 960, 1000, body)

def strip_board():
    rows = []
    for P in PALS:
        w, n = P["washi"], P["night"]
        def nav(p):
            return f'''<div style="display: flex; align-items: center; gap: 10px; height: 52px; padding: 0 16px; background: {p["ground"]}; border: 1px solid {p["line"]}; border-radius: 8px; width: 330px; box-sizing: border-box">
<div style="width: 30px; height: 30px">{F.svg(cols(p), viewBox=FVB, size=30, label="nav mark", prec=1)}</div>
<span style="{DISP}; font-size: 18px; font-weight: 800; letter-spacing: 0.05em; color: {p["ink"]}">GIDEON</span>
<a href="#" style="margin-left: auto; font-size: 13px; font-weight: 500; color: {p["link"]}; text-decoration: none">Guide</a><span style="font-size: 13px; color: {p["muted"]}">API</span>
<span style="display: inline-block; width: 30px; height: 30px; border-radius: 4px; background: {p["acc"]}"></span></div>'''
        rows.append(f'''<div style="display: flex; align-items: center; gap: 20px">
<div style="width: 150px; display: flex; flex-direction: column"><span style="{DISP}; font-size: 17px; font-weight: 800; color: #1d1b18">{P["name"]}</span><span style="font-size: 12px; color: #5f5a52">{P["kanji"]}</span></div>
{favs(w)}{nav(w)}<div style="background: {n["ground"]}; padding: 8px; border-radius: 10px; display: flex; gap: 14px; align-items: center">{favs(n)}{nav(n)}</div></div>''')
    body = f'''<div style="width: 1320px; height: 1000px; box-sizing: border-box; padding: 32px; background: #f3ede1; font-family: Inter, ui-sans-serif, system-ui, sans-serif; display: flex; flex-direction: column; gap: 18px">
<span style="{DISP}; font-size: 30px; line-height: 34px; font-weight: 800; letter-spacing: 0.04em; color: #1d1b18">Eight inks at real sizes</span>
{"".join(rows)}
</div>'''
    return page("GIDEON · Eight inks at real sizes", 1320, 1000, body)

if __name__ == "__main__":
    out = {f"Ink-{P['name'].replace(' ', '-')}.dc.html": ink_board(i, P) for i, P in enumerate(PALS)}
    out["Ground.dc.html"] = ground_board()
    out["Inks-Strip.dc.html"] = strip_board()
    for fn, t in out.items():
        open(os.path.join(ROOT, fn), "w").write(t); print(fn, len(t)//1024, "KB")
    for P in PALS:
        print(P["name"], {m: (P[m]["acc"], P[m]["cr"], P[m]["crink"], P[m]["crmut"]) for m in ("washi", "night")})
    json.dump(list(out), open("/tmp/final_files.json", "w"))
