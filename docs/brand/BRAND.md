# GIDEON · Ensō to sessen · Kogane 黄金

An ensō in ink, touched once by its tangent in gold; the carved G seal in gold, lower right; the Sumi family AMR
ground underneath, sumi-dashed (every cell edge its own brush dash, nodes left bare), in three levels refined around
the point of contact. GIDEON's own identity in the Sumi family (SETH, ADAM, FLUME): the same washi light ground, a
smoked-bronze dark ground in place of the family's silver ink, and gold in place of SETH's moegi.

This directory is the source of record. It is not part of the site (`srcExclude: brand/**`).

| File | What |
|---|---|
| `tokens.json` | every colour with its contrast on its ground, the fonts, the AMR levels |
| `generator/` | the Python that draws every stroke: `bundle.py` writes the marks, favicons and icons |

## Tokens

| Token | Light · washi | Dark · smoked bronze |
|---|---|---|
| ground | `#f3ede1` | `#11100c` |
| ink | `#1d1b18` (14.7:1) | `#efe9da` (15.7:1) |
| muted | `#5f5a52` (5.9:1) | `#bab29e` (9.0:1) |
| mesh | `#6f6a62` | `#948c79` |
| accent, gold (fills: seal, tangent, buttons) | `#855f00` (5.0:1) | `#f2c14e` (11.3:1) |
| links | `#765400` (5.9:1) | `#f3c861` (12.0:1) |
| text on the gold button | `#ffffff` (5.8:1) | `#1a1204` (11.1:1) |

Fonts: Shippori Mincho 800 (the name), Inter (text), JetBrains Mono (code).

## Where it is used

| Asset (`docs/public/`) | Used by |
|---|---|
| `logo.svg`, `logo-dark.svg` | the landing hero (`docs/index.md`), the README |
| `logo-nav.svg`, `logo-nav-dark.svg` | the navbar (`config.mts`, `themeConfig.logo`) |
| `favicon.svg`, `favicon-dark.svg`, `favicon-16.png`, `favicon-32.png`, `apple-touch-icon.png` | `config.mts`, `head` |
| `icons/stencils`, `icons/derivatives`, `icons/exact`, `icons/faces`, `icons/verified`, `icons/kernels` (each with `-dark`) | the six feature cards of the landing page |
| `icons/averages`, `icons/nabla`, `icons/licence` (each with `-dark`) | not used yet: drawn for the guide (averages and point values, the 3-D operators, the licence) |

The theme tokens are in `docs/.vitepress/theme/custom.css` (identity block first, layout after). The figures and the
casts of the documentation are drawn on the same two grounds: `docs/scripts/figure_theme.py` and
`docs/scripts/gideon_palette.py` carry the dark tokens above.

## Regenerate

```bash
cd docs/brand/generator
python3 bundle.py        # needs numpy, and rsvg-convert (librsvg) for the PNG favicons
```

It writes `out/gideon-identity/` (ignored; `GIDEON_BRAND_OUT` moves it): copy what changed into `docs/public/`.
Rerun as shipped, it reproduces every SVG of `docs/public` byte for byte. The only edit to the designer's generator:
its output paths, which pointed at the design session's machine, are now relative to this directory.
