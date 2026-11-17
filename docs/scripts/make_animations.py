#!/usr/bin/env python3
"""Render the animations of GIDEON's verification documentation from the frames of the verification programs.

Modelled on SETH's docs/scripts/make_doc_animations.py: each animation is drawn twice, on the two grounds of the docs
identity, washi (light) and smoked bronze (dark), with no title in the frame, a small time label under the panel, and no
margin of another colour. Per animation and theme it writes, into --out:

    <name>-<theme>.webm   VP9, constant quality (crf 36), 720 px wide
    <name>-<theme>.mp4    H.264, yuv420p, faststart
    <name>-<theme>.png    the poster: the last frame

Animations and the frames they read (written by src/tests/gideon_vv_advection.F90 in report mode, scripts/run_report.sh,
into <report>/r64): every curve is a solution advanced with GIDEON's compiled kernels.

    advection   advection-frames-{fv_d1_c2,fv_d1_c4,fv_d1_ul2}.dat   a wave packet over two periods, u_t + u_x = 0

Usage: make_animations.py [--report DIR] [--out DIR] [--only NAME ...] [--fps N] [--width PX]
Needs numpy, matplotlib and ffmpeg on PATH.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gideon_palette as P  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Theme:
    """The colours of one ground of the docs identity, and the line colours (gideon_palette)."""

    name: str
    bg: str
    bg_alt: str
    text: str
    mesh: str
    exact: str
    centred2: str
    centred4: str
    upwind: str


THEMES = (
    Theme('light', P.WASHI, P.WASHI_ALT, P.INK_LIGHT, P.MESH_LIGHT, P.INK_LIGHT, P.RIGHT, P.CENTRED, P.UPWIND),
    Theme('dark', P.BRONZE, P.BRONZE_ALT, P.INK_DARK, P.MESH_DARK, P.INK_DARK, P.RIGHT_D, P.CENTRED_D, P.UPWIND_D),
)
# kernel, label, the theme's colour attribute (the same roles as the static figure, make_figures.fig_advection)
ADVECTION = (('fv_d1_c2', 'centred c2 (fv_d1_c2)', 'centred2'), ('fv_d1_c4', 'centred c4 (fv_d1_c4)', 'centred4'),
             ('fv_d1_ul2', 'upwind ul2 (fv_d1_ul2)', 'upwind'))


def read_frames(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """x, t, u per frame, the exact solution per frame."""
    rows = [np.array(ln.split(), float) for ln in path.read_text().splitlines() if ln.strip()]
    x = rows[0]
    n = x.size
    data = np.array(rows[1:])
    return x, data[:, 0], data[:, 1:1 + n], data[:, 1 + n:1 + 2 * n]


def render(name: str, x: np.ndarray, t: np.ndarray, exact: np.ndarray,
           curves: list[tuple[str, str, np.ndarray]], ylabel: str, ylim: tuple[float, float], out: Path, fps: int,
           width: int) -> None:
    """One animation in both themes: frames to PNG, then ffmpeg to webm and mp4, and the poster."""
    for th in THEMES:
        tmp = Path(tempfile.mkdtemp(prefix=f'gideon-{name}-'))
        try:
            with mpl.rc_context({'figure.facecolor': th.bg, 'savefig.facecolor': th.bg, 'axes.facecolor': th.bg_alt,
                                 'axes.edgecolor': th.mesh, 'axes.labelcolor': th.text, 'text.color': th.text,
                                 'xtick.color': th.mesh, 'ytick.color': th.mesh, 'xtick.labelcolor': th.text,
                                 'ytick.labelcolor': th.text, 'grid.color': th.mesh, 'legend.facecolor': th.bg,
                                 'legend.edgecolor': th.mesh, 'font.size': 9}):
                for f in range(t.size):
                    fig = plt.figure(figsize=(7.2, 3.2), dpi=100)
                    ax = fig.add_axes([0.1, 0.2, 0.88, 0.76])
                    ax.plot(x, exact[f], '-', color=th.exact, lw=1.0, label='exact')
                    for lab, attr, field in curves:
                        ax.plot(x, field[f], '-', color=getattr(th, attr), lw=1.2, label=lab)
                    ax.set_xlim(0.0, 1.0)
                    ax.set_ylim(*ylim)
                    ax.set_ylabel(ylabel)
                    ax.grid(alpha=0.25)
                    ax.legend(loc='upper right', fontsize=7.5, frameon=True)
                    fig.text(0.1, 0.04, f't = {t[f]:.3f} (periods)', fontsize=9, color=th.text)
                    fig.savefig(tmp / f'f{f:04d}.png')
                    plt.close(fig)
            out.mkdir(parents=True, exist_ok=True)
            base = out / f'{name}-{th.name}'
            scale = f'scale={width}:-2:flags=lanczos'
            subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(tmp / 'f%04d.png'),
                            '-vf', scale, '-c:v', 'libvpx-vp9', '-crf', '36', '-b:v', '0', '-pix_fmt', 'yuv420p',
                            f'{base}.webm'], check=True)
            subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', str(tmp / 'f%04d.png'),
                            '-vf', scale, '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-movflags', '+faststart',
                            f'{base}.mp4'], check=True)
            subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-i', str(tmp / f'f{t.size - 1:04d}.png'), '-vf',
                            scale, f'{base}.png'], check=True)
            print(f'  {base.relative_to(ROOT) if base.is_relative_to(ROOT) else base}.{{webm,mp4,png}}')
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--report', type=Path, default=ROOT / 'exe' / 'report')
    ap.add_argument('--out', type=Path, default=ROOT / 'docs' / 'public' / 'verification' / 'casts')
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--fps', type=int, default=15)
    ap.add_argument('--width', type=int, default=720)
    args = ap.parse_args()
    if shutil.which('ffmpeg') is None:
        print('ffmpeg not found: no animation rendered')
        return 1
    if not args.only or 'advection' in args.only:
        files = [args.report / 'r64' / f'advection-frames-{k}.dat' for k, _, _ in ADVECTION]
        if not all(p.exists() for p in files):
            print('advection: missing frames (scripts/run_report.sh writes them)')
            return 1
        curves, x, t, exact = [], None, None, None
        for (_, lab, attr), p in zip(ADVECTION, files, strict=True):
            x, t, u, ex = read_frames(p)
            curves.append((lab, attr, u))
            exact = ex
        lim = 1.15 * max(np.abs(exact).max(), max(np.abs(c[2]).max() for c in curves))
        print('advection')
        render('advection', x, t, exact, curves, 'u (cell averages)', (-lim, lim), args.out, args.fps, args.width)
    return 0


if __name__ == '__main__':
    sys.exit(main())
