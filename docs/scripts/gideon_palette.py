"""The colours of GIDEON's documentation figures and animations: the Sumi family identity applied to data.

Modelled on SETH's docs/scripts/seth_palette.py. GIDEON wears SETH's Neon Moegi ink until it has its own identity
(docs/.vitepress/theme/custom.css), so the grounds, the inks and the four categorical data colours are SETH's, which
were validated there with the dataviz skill's validate_palette.js (lightness band, chroma floor, colour-blind and
normal-vision separation of adjacent pairs, contrast against washi #f3ede1 and silver ink #0e1116). The colours live
here once, for the figures (make_figures.py) and the animations (make_animations.py).

The static figures are drawn twice, with the light colours and, inside figure_theme.theme('dark') (which maps the
neutral greys and the grounds), with the dark ones; the animations draw each theme with its own colours.
"""

from __future__ import annotations

# identity grounds (custom.css)
WASHI, WASHI_ALT, INK_LIGHT, MESH_LIGHT = '#f3ede1', '#ebe3d3', '#1d1b18', '#6f6a62'
SILVER, SILVER_ALT, INK_DARK, MESH_DARK = '#0e1116', '#0a0d11', '#e8ecf0', '#8b949e'

# categorical data colours, in fixed order (never cycled): blue, moegi, plum, umber
BLUE, MOEGI, PLUM, UMBER = '#1f63a8', '#4a7a06', '#8a4f9e', '#a8661c'
BLUE_D, MOEGI_D, PLUM_D, UMBER_D = '#3f8fd6', '#6aa41c', '#b06ac4', '#c27a2a'
LINES_LIGHT = (BLUE, MOEGI, PLUM, UMBER)
LINES_DARK = (BLUE_D, MOEGI_D, PLUM_D, UMBER_D)

# the roles of the figures: centred kernels blue, upwind (left side) moegi, the right side plum, ADAM umber
CENTRED, UPWIND, RIGHT, ADAM = BLUE, MOEGI, PLUM, UMBER
CENTRED_D, UPWIND_D, RIGHT_D, ADAM_D = BLUE_D, MOEGI_D, PLUM_D, UMBER_D

# a sequential ramp per family, light to dark with the order (the order is the shade, the family the hue)
RAMP_BLUE = ('#8db6e0', '#4f8fcf', '#1f63a8', '#174d84', '#0f3660')
RAMP_MOEGI = ('#a5c76a', '#76a52a', '#4a7a06', '#3a6004', '#2a4603')
RAMP_BLUE_D = ('#2a5f93', '#3f8fd6', '#6aa8e6', '#9cc6f0', '#c8dff6')
RAMP_MOEGI_D = ('#3e6a10', '#6aa41c', '#8fc23e', '#b5d77a', '#d3e8ad')


def ramp(base: tuple[str, ...], k: int, n: int) -> str:
    """The k-th of n shades of a ramp (k = 0 the lightest on washi), spread over the ramp."""
    if n <= 1:
        return base[len(base) // 2]
    return base[round(k * (len(base) - 1) / (n - 1))]
