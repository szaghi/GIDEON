import json, os, html
from final import PALS, M, F, FVB, ICN, cols, mark, sw, page, DISP, MONO, TAG
from icons import ICONS
from pal import cr

ROOT = os.path.join(os.environ.get("GIDEON_BRAND_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")), "boards")  # was /home/claude/gideon/project
K = [p for p in PALS if p["key"] == "kogane"][0]
W, N = dict(K["washi"]), dict(K["night"])
W.update(pressed="#6b4c00", on="#ffffff", focus="#1d1b18")
N.update(pressed="#c99a2e", on="#1a1204", focus="#f3c861")

TOKENS = [
 ("ground", "Page background: washi paper, smoked bronze at night."),
 ("alt", "Sidebars, code blocks, alternating bands."),
 ("raised", "Cards, menus, the favicon tile."),
 ("line", "Hairlines and dividers, decorative only."),
 ("ink", "Body text and the ensō."),
 ("muted", "Secondary text and captions."),
 ("mesh", "The sumi-dashed AMR ground; control borders."),
 ("acc", "Kogane: the tangent, the seal, primary button fill."),
 ("link", "Links and the active nav item."),
 ("pressed", "Primary button hover and pressed."),
 ("on", "Text and icons on the gold fill."),
 ("focus", "2px keyboard focus ring with a 2px ground gap."),
]
NAMES = dict(ground="ground", alt="ground-alt", raised="ground-raised", line="line", ink="ink", muted="ink-muted", mesh="mesh",
             acc="seal / accent", link="link", pressed="action-pressed", on="on-action", focus="focus-ring")

def ratio(p, k):
    if k in ("ground", "alt", "raised", "line"): return ""
    if k == "on": return f"{cr(p['on'], p['acc']):.1f}:1 on accent"
    return f"{cr(p[k], p['ground']):.1f}:1"

def identity():
    def col(p, mode):
        ladder = "".join(
            f'<div style="display: flex; flex-direction: column; align-items: center; gap: 6px"><div style="width: {px}px; height: {px}px; border-radius: {max(2, round(px * 0.19))}px; background: {p["raised"]}; border: 1px solid {p["line"]}; display: flex; align-items: center; justify-content: center">{F.svg(cols(p), viewBox=FVB, size=round(px * 0.86), label=f"favicon {px}px", prec=1)}</div><span style="{MONO}; font-size: 10px; color: {p["muted"]}">{px}</span></div>'
            for px in (96, 64, 48, 32, 24, 16))
        rows = "".join(
            f'<div style="display: grid; grid-template-columns: 22px 120px 70px 1fr; gap: 8px; align-items: center"><span style="width: 22px; height: 16px; border-radius: 4px; background: {p[k]}; border: 1px solid {p["line"]}"></span><span style="{MONO}; font-size: 11px; color: {p["ink"]}">{NAMES[k]}</span><span style="{MONO}; font-size: 11px; color: {p["muted"]}">{p[k]}</span><span style="font-size: 11px; color: {p["muted"]}">{ratio(p, k)}</span></div>'
            for k, _ in TOKENS)
        icons = "".join(f'<div style="display: flex; flex-direction: column; gap: 5px; align-items: center; width: 44px">{ICN[k].svg(cols(p), viewBox="0 0 24 24", size=30, label=ICONS[k], prec=1)}</div>' for k in range(9))
        return f'''<div style="width: 480px; box-sizing: border-box; padding: 28px 32px; background: {p["ground"]}; color: {p["ink"]}; display: flex; flex-direction: column; gap: 18px">
<div style="display: flex; justify-content: space-between; align-items: baseline"><p style="margin: 0; font-size: 12px; font-weight: 500; letter-spacing: 0.02em; color: {p["muted"]}">{mode}</p>{sw(p)}</div>
<div style="width: 416px; height: 416px">{mark(p, 416, "GIDEON logo, Ensō to sessen, Kogane, " + mode)}</div>
<div style="display: flex; flex-direction: column; gap: 6px"><span style="{DISP}; font-size: 44px; line-height: 44px; font-weight: 800; letter-spacing: 0.04em">GIDEON</span><span style="font-size: 12px; line-height: 16px; font-weight: 500; letter-spacing: 0.02em; color: {p["muted"]}">{TAG}</span></div>
<div style="display: flex; gap: 12px; align-items: flex-end">{ladder}</div>
<div style="display: flex; justify-content: space-between">{icons}</div>
<div style="display: flex; flex-direction: column; gap: 6px">{rows}</div>
</div>'''
    body = f'''<div style="width: 960px; height: 1240px; display: flex; font-family: Inter, ui-sans-serif, system-ui, sans-serif">{col(W, "Washi")}{col(N, "Kogane night")}</div>'''
    return page("GIDEON · Kogane identity", 960, 1240, body)

FEATURES = [
 (0, "Exact coefficients", "Every coefficient is derived in exact rational arithmetic and cross-checked against ADAM's tables. Nothing is typed by hand."),
 (4, "243 kernels", "Derivatives d1 to d6, face values, face derivatives and cell averages, with centred or upwind stencils of 1 to 9 cells."),
 (5, "Callable from inside your loops", "One pure, device-callable kernel per operator and order, at 0.17 to 0.82 of ADAM's cost on a GPU and a CPU core."),
 (6, "Higher-level helpers", "Run-time order dispatchers and 3-D gradient, divergence, curl and Laplacian operators."),
 (7, "Verified", "Exactness checks, observed orders in quadruple precision, bitwise equivalence with ADAM, 22 of 22 mutants killed."),
 (8, "Multi-licensed", "Choose GPL-3.0-or-later, BSD-2-Clause, BSD-3-Clause or MIT."),
]

def docs(p, mode):
    feats = "".join(f'''<div style="background: {p["raised"]}; border: 1px solid {p["line"]}; border-radius: 8px; padding: 20px; display: flex; flex-direction: column; gap: 10px">
<div style="width: 44px; height: 44px; border-radius: 8px; background: {p["alt"]}; display: flex; align-items: center; justify-content: center">{ICN[i].svg(cols(p), viewBox="0 0 24 24", size=30, label="", prec=1)}</div>
<h3 style="margin: 0; font-size: 16px; line-height: 22px; font-weight: 600">{t}</h3>
<p style="margin: 0; font-size: 14px; line-height: 22px; color: {p["muted"]}">{d}</p></div>''' for i, t, d in FEATURES)
    nav = "".join(f'<a href="#" style="font-size: 14px; font-weight: 500; color: {p["link"] if k == 0 else p["ink"]}; text-decoration: none">{t}</a>' for k, t in enumerate(("Guide", "Tutorials", "Cookbook", "Verification", "API", "GitHub")))
    body = f'''<div style="min-height: 900px; background: {p["ground"]}; color: {p["ink"]}; font-family: Inter, ui-sans-serif, system-ui, sans-serif">
<header style="border-bottom: 1px solid {p["line"]}"><div style="max-width: 1152px; margin: 0 auto; padding: 0 24px; height: 64px; display: flex; align-items: center; gap: 12px; flex-wrap: wrap">
<div style="width: 32px; height: 32px">{F.svg(cols(p), viewBox=FVB, size=32, label="GIDEON", prec=1)}</div>
<span style="{DISP}; font-size: 20px; font-weight: 800; letter-spacing: 0.05em">GIDEON</span>
<nav style="margin-left: auto; display: flex; gap: 22px; flex-wrap: wrap">{nav}</nav></div></header>
<main style="max-width: 1152px; margin: 0 auto; padding: 56px 24px 64px; display: flex; flex-direction: column; gap: 56px">
<section style="display: flex; gap: 48px; align-items: center; flex-wrap: wrap">
<div style="flex: 1 1 420px; display: flex; flex-direction: column; gap: 18px">
<h1 style="margin: 0; {DISP}; font-size: 56px; line-height: 60px; font-weight: 800; letter-spacing: 0.04em; color: {p["ink"]}">GIDEON</h1>
<p style="margin: 0; font-size: 22px; line-height: 30px; font-weight: 600; color: {p["link"]}">{TAG}</p>
<p style="margin: 0; font-size: 17px; line-height: 27px; color: {p["muted"]}; max-width: 520px">Finite-difference and finite-volume stencil operators for high-performance computing, in modern Fortran, on CPUs and GPUs.</p>
<div style="display: flex; gap: 12px; flex-wrap: wrap; padding-top: 6px">
<a href="#" style="display: inline-flex; align-items: center; height: 44px; padding: 0 22px; border-radius: 4px; background: {p["acc"]}; color: {p["on"]}; font-size: 15px; font-weight: 600; text-decoration: none">Get started</a>
<a href="#" style="display: inline-flex; align-items: center; height: 44px; padding: 0 22px; border-radius: 4px; border: 1px solid {p["mesh"]}; color: {p["ink"]}; font-size: 15px; font-weight: 600; text-decoration: none">Tutorials</a>
<a href="#" style="display: inline-flex; align-items: center; height: 44px; padding: 0 22px; border-radius: 4px; border: 1px solid {p["mesh"]}; color: {p["ink"]}; font-size: 15px; font-weight: 600; text-decoration: none">View on GitHub</a></div>
</div>
<div style="flex: 0 1 360px; display: flex; justify-content: center">{mark(p, 340, "GIDEON logo")}</div>
</section>
<section style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 16px">{feats}</section>
</main>
<footer style="border-top: 1px solid {p["line"]}"><div style="max-width: 1152px; margin: 0 auto; padding: 20px 24px; font-size: 13px; color: {p["muted"]}">GIDEON · GPL-3.0-or-later, BSD-2-Clause, BSD-3-Clause or MIT</div></footer>
</div>'''
    return page(f"GIDEON docs · {mode}", 1280, 1100, body, bg=p["ground"], link=p["link"])

CSS = f"""/* .vitepress/theme/custom.css — GIDEON · Ensō to sessen · Kogane 黄金 */
@import url('https://fonts.googleapis.com/css2?family=Shippori+Mincho:wght@800&family=Inter:wght@400;500;600&family=JetBrains+Mono&display=swap');

:root {{
  --vp-font-family-base: Inter, ui-sans-serif, system-ui, sans-serif;
  --vp-font-family-mono: 'JetBrains Mono', ui-monospace, monospace;

  --vp-c-bg: {W['ground']};
  --vp-c-bg-alt: {W['alt']};
  --vp-c-bg-elv: {W['raised']};
  --vp-c-bg-soft: {W['alt']};
  --vp-c-divider: {W['line']};
  --vp-c-border: {W['mesh']};
  --vp-c-text-1: {W['ink']};
  --vp-c-text-2: {W['muted']};

  --vp-c-brand-1: {W['link']};   /* links, active nav: {cr(W['link'], W['ground']):.1f}:1 */
  --vp-c-brand-2: {W['acc']};
  --vp-c-brand-3: {W['acc']};    /* button fill */
  --vp-c-brand-soft: rgba(133, 95, 0, 0.14);

  --vp-button-brand-bg: {W['acc']};
  --vp-button-brand-text: {W['on']};
  --vp-button-brand-hover-bg: {W['pressed']};
  --vp-button-brand-active-bg: {W['pressed']};

  --vp-home-hero-name-color: {W['ink']};
  --gideon-mesh: {W['mesh']};
  --gideon-focus: {W['focus']};
}}

.dark {{
  --vp-c-bg: {N['ground']};
  --vp-c-bg-alt: {N['alt']};
  --vp-c-bg-elv: {N['raised']};
  --vp-c-bg-soft: {N['raised']};
  --vp-c-divider: {N['line']};
  --vp-c-border: {N['mesh']};
  --vp-c-text-1: {N['ink']};
  --vp-c-text-2: {N['muted']};

  --vp-c-brand-1: {N['link']};
  --vp-c-brand-2: {N['acc']};
  --vp-c-brand-3: {N['acc']};
  --vp-c-brand-soft: rgba(242, 193, 78, 0.16);

  --vp-button-brand-bg: {N['acc']};
  --vp-button-brand-text: {N['on']};
  --vp-button-brand-hover-bg: {N['pressed']};
  --vp-button-brand-active-bg: {N['pressed']};

  --vp-home-hero-name-color: {N['ink']};
  --gideon-mesh: {N['mesh']};
  --gideon-focus: {N['focus']};
}}

.VPHero .name, .VPNavBarTitle .title {{
  font-family: 'Shippori Mincho', 'Hiragino Mincho ProN', serif;
  font-weight: 800;
  letter-spacing: 0.04em;
}}

:focus-visible {{
  outline: 2px solid var(--gideon-focus);
  outline-offset: 2px;
}}

/* config.mts → head: ['meta', {{ name: 'theme-color', content: '{W['ground']}', media: '(prefers-color-scheme: light)' }}],
                       ['meta', {{ name: 'theme-color', content: '{N['ground']}', media: '(prefers-color-scheme: dark)' }}] */
"""

def vp():
    body = f'''<div style="width: 760px; height: 1240px; box-sizing: border-box; padding: 28px 32px; background: {W["alt"]}; font-family: Inter, ui-sans-serif, system-ui, sans-serif; display: flex; flex-direction: column; gap: 14px">
<span style="{DISP}; font-size: 26px; line-height: 30px; font-weight: 800; letter-spacing: 0.04em; color: {W["ink"]}">VitePress override</span>
<span style="font-size: 13px; line-height: 19px; color: {W["muted"]}">Drop into .vitepress/theme/custom.css. The night ground moves from #0e1116 to Kogane's smoked bronze {N["ground"]}, so update the theme-color meta too.</span>
<pre style="margin: 0; padding: 16px 18px; background: {W["raised"]}; border: 1px solid {W["line"]}; border-radius: 8px; {MONO}; font-size: 11.5px; line-height: 17px; color: {W["ink"]}; white-space: pre-wrap; overflow: hidden">{html.escape(CSS)}</pre>
</div>'''
    return page("GIDEON · VitePress override", 760, 1240, body, bg=W["alt"])

if __name__ == "__main__":
    out = {"Identity-Kogane.dc.html": identity(), "Docs-Washi.dc.html": docs(W, "Washi"), "Docs-Night.dc.html": docs(N, "Night"), "VitePress.dc.html": vp()}
    for fn, t in out.items():
        open(os.path.join(ROOT, fn), "w").write(t); print(fn, len(t) // 1024, "KB")
    open(os.path.join(ROOT, "custom.css"), "w").write(CSS)  # was /home/claude/gideon/custom.css
