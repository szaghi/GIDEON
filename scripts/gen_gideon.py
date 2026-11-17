#!/usr/bin/env python3
"""GIDEON generator: every stencil row derived in exact rational arithmetic by its meaning, then emitted as Fortran.

A row is the set of weights w_k such that sum_k w_k L_k(p) = T(p) for every polynomial p of degree < N (N cells or
nodes in the window): L_k samples the data (a point value, or a cell average), T is the target (a point value or
derivative, or the cell average of a derivative). The system is solved exactly (fractions.Fraction); the row's
maximum exact degree K is then measured on the monomials beyond N, and its order of accuracy is K + 1 - n for the
n-th derivative. Nothing is copied from a table: ADAM's committed tables (notes/adam-reference/) are parsed only to
cross-check the derived rows, and a mismatch is reported, never silently adopted.

Geometry, in units of the grid step (h = 1); windows are offsets, see docs/guide/conventions.md:

    fd node : point values at x = o, the target at the node o = 0
    fd face : point values at x = o - 1/2, the target at the face x = 0 between nodes 0 and 1
    fv face : cell o spans [o - 1, o], the target at the face x = 0 between cells 0 and 1 (SETH's face)
    fv cell : cell o spans [o - 1/2, o + 1/2], the target over cell 0

Left-upwind rows are never stored: they are the right-upwind rows mirrored (times (-1)^n), which the generator
proves by solving them independently and comparing.

Usage:
    python3 scripts/gen_gideon.py --check     derive, verify and cross-check; write nothing
    python3 scripts/gen_gideon.py             also (re)write the generated Fortran sources
    python3 scripts/gen_gideon.py --check --json FILE
                                              also write the derived rows, windows, orders and findings as JSON (the
                                              data of the documentation's figures and tables, docs/scripts)
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from decimal import Decimal, getcontext
from fractions import Fraction as Fr
from math import factorial
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ADAM_REF = ROOT / 'notes' / 'adam-reference' / 'adam_fdv_operators_library.F90'

S_CENTRED = range(1, 6)  # half-stencils of the centred families (ADAM's S_MAX = 5), owner decision D2
S_UPWIND = range(1, 10)  # stencil widths of the upwind families, D2
NDER_MAX = 6  # highest derivative


# ---------------------------------------------------------------------------------------------------------------------
# exact algebra
# ---------------------------------------------------------------------------------------------------------------------
def solve(a: list[list[Fr]], b: list[Fr]) -> list[Fr]:
    """Solve a x = b exactly (Gauss-Jordan with pivoting on non-zero entries)."""
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = next((r for r in range(c, n) if m[r][c] != 0), None)
        if p is None:
            raise ValueError('singular system')
        m[c], m[p] = m[p], m[c]
        piv = m[c][c]
        m[c] = [v / piv for v in m[c]]
        for r in range(n):
            if r != c and m[r][c] != 0:
                f = m[r][c]
                m[r] = [vr - f * vc for vr, vc in zip(m[r], m[c], strict=True)]
    return [m[r][n] for r in range(n)]


def mono_point(x: Fr, j: int) -> Fr:
    """Point value of x^j."""
    return x**j


def mono_avg(a: Fr, b: Fr, j: int) -> Fr:
    """Average of x^j over [a, b]."""
    return (b ** (j + 1) - a ** (j + 1)) / ((j + 1) * (b - a))


def mono_deriv_point(n: int, x: Fr, j: int) -> Fr:
    """n-th derivative of x^j at x."""
    if j < n:
        return Fr(0)
    return Fr(factorial(j), factorial(j - n)) * x ** (j - n)


def mono_deriv_avg(n: int, a: Fr, b: Fr, j: int) -> Fr:
    """Average over [a, b] of the n-th derivative of x^j."""
    if j < n:
        return Fr(0)
    return Fr(factorial(j), factorial(j - n)) * mono_avg(a, b, j - n)


# ---------------------------------------------------------------------------------------------------------------------
# operators: sampling (data) and target functionals
# ---------------------------------------------------------------------------------------------------------------------
HALF = Fr(1, 2)


def sampler(data: str, place: str, o: int):
    """Functional L_o of the datum at offset o, as a function of the monomial degree j."""
    if data == 'fd':
        x = Fr(o) if place in ('node', 'cell') else Fr(o) - HALF
        return lambda j: mono_point(x, j)
    if place == 'face':  # cell o spans [o - 1, o]
        a, b = Fr(o - 1), Fr(o)
    else:  # cell o spans [o - 1/2, o + 1/2]
        a, b = Fr(o) - HALF, Fr(o) + HALF
    return lambda j: mono_avg(a, b, j)


def target(kind: str, n: int):
    """Target functional T as a function of the monomial degree j."""
    if kind == 'point':  # n-th derivative at x = 0 (node or face)
        return lambda j: mono_deriv_point(n, Fr(0), j)
    if kind == 'cellavg':  # average over cell 0 of the n-th derivative
        return lambda j: mono_deriv_avg(n, -HALF, HALF, j)
    raise ValueError(kind)


@dataclass
class Row:
    """A derived row: weights by offset, its window, its exactness."""

    weights: dict[int, Fr]
    lo: int
    hi: int
    n: int
    exact_degree: int

    @property
    def order(self) -> int:
        return self.exact_degree + 1 - self.n

    def is_zero(self) -> bool:
        return all(w == 0 for w in self.weights.values())


def exact_degree(apply, tgt, jmax: int) -> int:
    """Largest K such that the operator is exact on every monomial of degree <= K (checked up to jmax)."""
    for j in range(jmax + 1):
        if apply(j) != tgt(j):
            return j - 1
    return jmax


def derive(data: str, place: str, kind: str, n: int, offsets: list[int]) -> Row:
    """Derive the row of the operator on the given window by its meaning."""
    samp = [sampler(data, place, o) for o in offsets]
    tgt = target(kind, n)
    npts = len(offsets)
    a = [[samp[k](j) for k in range(npts)] for j in range(npts)]
    b = [tgt(j) for j in range(npts)]
    w = solve(a, b)
    apply = lambda j: sum(wk * s(j) for wk, s in zip(w, samp, strict=True))  # noqa: E731
    k = exact_degree(apply, tgt, npts + 4)
    if k < npts - 1:
        raise AssertionError(f'{data} {place} {kind} n={n} {offsets}: not exact on its own system')
    return Row(dict(zip(offsets, w, strict=True)), min(offsets), max(offsets), n, k)


def composite_order(data: str, n: int, weights: dict[int, Fr]) -> int:
    """Order of a cell-average derivative given by its full weights on cells (fv cell geometry)."""
    offs = sorted(weights)
    samp = {o: sampler(data, 'cell', o) for o in offs}
    tgt = target('cellavg', n)
    apply = lambda j: sum(weights[o] * samp[o](j) for o in offs)  # noqa: E731
    return exact_degree(apply, tgt, len(offs) + 6) + 1 - n


# ---------------------------------------------------------------------------------------------------------------------
# the families
# ---------------------------------------------------------------------------------------------------------------------
@dataclass
class Kernel:
    """One generated kernel: name, family, its rows and window, its order."""

    name: str
    family: str  # fd_face_c, fd_node_c, fd_facedn_c, fd_face_u, fd_face_b, fv_face_c, fv_facedn_c, fv_cell_c,
    #               fv_face_u, fv_facedn_u, fv_cell_u, deconv
    data: str
    op: str  # face, dN, face_dN, point, average
    bias: str  # c, ur, ul, br, bl
    n: int
    s: int  # half-stencil (centred), width (upwind) or k (biased: k+1 points on one side, k on the other)
    order: int
    lo: int
    hi: int
    parity: int = 0  # centred: +1 symmetric pairs, -1 antisymmetric pairs
    row: dict[int, Fr] = field(default_factory=dict)  # stored row (centred: pair half-row incl. c0; upwind: ur row;
    #                                                     biased: br row)
    sumabs: Fr = Fr(0)  # sum of |weights| of the full operator (round-off scale of the tests)
    full: dict[int, Fr] = field(default_factory=dict)  # every weight of the operator by offset (orders check)


def check_pairs(row: Row, centre: str) -> tuple[int, dict[int, Fr]]:
    """Pair form of a centred row: parity and the half-row (key 0 = centre weight for node/cell rows)."""
    w = row.weights
    if centre == 'node':
        s = row.hi
        sym = all(w[m] == w[-m] for m in range(1, s + 1))
        anti = all(w[m] == -w[-m] for m in range(1, s + 1)) and w[0] == 0
        half = {m: w[m] for m in range(0 if sym and not anti else 1, s + 1)}
    else:  # face: pairs (m, 1 - m)
        s = row.hi
        sym = all(w[m] == w[1 - m] for m in range(1, s + 1))
        anti = all(w[m] == -w[1 - m] for m in range(1, s + 1))
        half = {m: w[m] for m in range(1, s + 1)}
    if sym == anti:
        raise AssertionError(f'row is neither symmetric nor antisymmetric: {w}')
    return (1 if sym else -1), half


def build() -> list[Kernel]:
    """Derive every kernel of P1."""
    ks: list[Kernel] = []

    def add(k: Kernel) -> None:
        if any(x.name == k.name for x in ks):
            raise AssertionError(f'duplicate kernel name {k.name}')
        ks.append(k)

    # fd centred: face interpolation and node derivatives
    for s in S_CENTRED:
        r = derive('fd', 'face', 'point', 0, list(range(1 - s, s + 1)))
        par, half = check_pairs(r, 'face')
        add(Kernel(f'gideon_fd_face_c{r.order}', 'fd_face_c', 'fd', 'face', 'c', 0, s, r.order, 1 - s, s, par, half,
                   sum(abs(v) for v in r.weights.values()), dict(r.weights)))
    for n in range(1, NDER_MAX + 1):
        for s in S_CENTRED:
            if 2 * s + 1 <= n:
                continue
            r = derive('fd', 'node', 'point', n, list(range(-s, s + 1)))
            if r.is_zero() or r.order < 1:
                continue
            par, half = check_pairs(r, 'node')
            add(Kernel(f'gideon_fd_d{n}_c{r.order}', 'fd_node_c', 'fd', f'd{n}', 'c', n, s, r.order, -s, s, par, half,
                       sum(abs(v) for v in r.weights.values()), dict(r.weights)))

    # fv centred: face value, face derivatives, cell-average derivatives (flux difference of the face rows)
    face_rows: dict[tuple[int, int], Row] = {}
    for n in range(0, NDER_MAX):
        for s in S_CENTRED:
            if 2 * s <= n:
                continue
            r = derive('fv', 'face', 'point', n, list(range(1 - s, s + 1)))
            if r.is_zero() or r.order < 1:
                continue
            face_rows[(n, s)] = r
            par, half = check_pairs(r, 'face')
            op = 'face' if n == 0 else f'face_d{n}'
            fam = 'fv_face_c' if n == 0 else 'fv_facedn_c'
            add(Kernel(f'gideon_fv_{op}_c{r.order}', fam, 'fv', op, 'c', n, s, r.order, 1 - s, s, par, half,
                       sum(abs(v) for v in r.weights.values()), dict(r.weights)))
    for n in range(1, NDER_MAX + 1):
        for s in S_CENTRED:
            if (n - 1, s) not in face_rows:
                continue
            g = face_rows[(n - 1, s)]
            full = {o: Fr(0) for o in range(-s, s + 1)}
            for o, w in g.weights.items():  # g_{+1/2} on cells o (face between 0 and 1), g_{-1/2} on cells o - 1
                full[o] += w
                full[o - 1] -= w
            order = composite_order('fv', n, full)
            if order < 1:
                continue
            par, half = check_pairs(g, 'face')
            add(Kernel(f'gideon_fv_d{n}_c{order}', 'fv_cell_c', 'fv', f'd{n}', 'c', n, s, order, -s, s, par, half,
                       sum(abs(v) for v in full.values()), full))

    # fv upwind: right-side stencils derived, left-side ones derived independently and proven to be their mirrors
    up_rows: dict[tuple[int, int], Row] = {}
    for n in range(0, NDER_MAX):
        for s in S_UPWIND:
            if s <= n:
                continue
            r = derive('fv', 'face', 'point', n, list(range(1, s + 1)))
            left = derive('fv', 'face', 'point', n, list(range(1 - s, 1)))
            sign = -1 if n % 2 else 1
            for o in range(1, s + 1):
                if left.weights[1 - o] != sign * r.weights[o]:
                    raise AssertionError(f'fv face n={n} s={s}: the left row is not the mirror of the right one')
            if r.is_zero() or r.order < 1:
                continue
            up_rows[(n, s)] = r
            op = 'face' if n == 0 else f'face_d{n}'
            fam = 'fv_face_u' if n == 0 else 'fv_facedn_u'
            for bias, lo, hi in (('ur', 1, s), ('ul', 1 - s, 0)):
                fw = dict(r.weights) if bias == 'ur' else {1 - o: sign * w for o, w in r.weights.items()}
                add(Kernel(f'gideon_fv_{op}_{bias}{r.order}', fam, 'fv', op, bias, n, s, r.order, lo, hi, 0,
                           dict(r.weights), sum(abs(v) for v in r.weights.values()), fw))
    for n in range(1, NDER_MAX + 1):
        for s in S_UPWIND:
            if (n - 1, s) not in up_rows:
                continue
            g = up_rows[(n - 1, s)]
            full = {o: Fr(0) for o in range(0, s + 1)}
            for o, w in g.weights.items():  # ur: face +1/2 on cells 1..s, face -1/2 on cells 0..s-1
                full[o] += w
                full[o - 1] -= w
            order = composite_order('fv', n, full)
            mirror = {-o: (-1) ** n * w for o, w in full.items()}
            if composite_order('fv', n, mirror) != order:
                raise AssertionError(f'fv cell d{n} s={s}: the mirrored left operator does not keep the order')
            if order < 1:
                continue
            for bias, lo, hi in (('ur', 0, s), ('ul', -s, 0)):
                fw = full if bias == 'ur' else mirror
                add(Kernel(f'gideon_fv_d{n}_{bias}{order}', 'fv_cell_u', 'fv', f'd{n}', bias, n, s, order, lo, hi, 0,
                           dict(g.weights), sum(abs(v) for v in full.values()), fw))

    # deconvolution: cell average -> point value at the centre, point values -> cell average
    for s in S_CENTRED:
        r = derive('fv', 'cell', 'point', 0, list(range(-s, s + 1)))
        par, half = check_pairs(r, 'node')
        add(Kernel(f'gideon_fv_point_c{r.order}', 'deconv', 'fv', 'point', 'c', 0, s, r.order, -s, s, par, half,
                   sum(abs(v) for v in r.weights.values()), dict(r.weights)))
        r = derive('fd', 'cell', 'cellavg', 0, list(range(-s, s + 1)))
        par, half = check_pairs(r, 'node')
        add(Kernel(f'gideon_fd_average_c{r.order}', 'deconv', 'fd', 'average', 'c', 0, s, r.order, -s, s, par, half,
                   sum(abs(v) for v in r.weights.values()), dict(r.weights)))

    # phase 2 (appended: the first release's kernels keep their order, and their golden values their lines)
    # 2.1 fd centred face derivatives from point values: FLUME's viscous rows, (27(W1-W0) - (W2-W-1))/24 at order 4
    for n in range(1, NDER_MAX):
        for s in S_CENTRED:
            if 2 * s <= n:
                continue
            r = derive('fd', 'face', 'point', n, list(range(1 - s, s + 1)))
            if r.is_zero() or r.order < 1:
                continue
            par, half = check_pairs(r, 'face')
            add(Kernel(f'gideon_fd_face_d{n}_c{r.order}', 'fd_facedn_c', 'fd', f'face_d{n}', 'c', n, s, r.order, 1 - s,
                       s, par, half, sum(abs(v) for v in r.weights.values()), dict(r.weights)))
    # 2.2 fd one-sided face values (extrapolation by half a step): right side derived, left side derived
    # independently and proven to be its mirror (points at x = o - 1/2: the mirror of offset o is 1 - o)
    for s in S_UPWIND:
        r = derive('fd', 'face', 'point', 0, list(range(1, s + 1)))
        left = derive('fd', 'face', 'point', 0, list(range(1 - s, 1)))
        for o in range(1, s + 1):
            if left.weights[1 - o] != r.weights[o]:
                raise AssertionError(f'fd face s={s}: the left row is not the mirror of the right one')
        for bias, lo, hi in (('ur', 1, s), ('ul', 1 - s, 0)):
            fw = dict(r.weights) if bias == 'ur' else {1 - o: w for o, w in r.weights.items()}
            add(Kernel(f'gideon_fd_face_{bias}{r.order}', 'fd_face_u', 'fd', 'face', bias, 0, s, r.order, lo, hi, 0,
                       dict(r.weights), sum(abs(v) for v in r.weights.values()), fw))
    # 2.2 fd upwind-biased face values (interpolation): k + 1 points on one side of the face, k on the other;
    # br (offsets 1-k..k+1, one more on the right) derived, bl (offsets -k..k) derived and proven its mirror
    for k in range(1, (max(S_UPWIND) - 1) // 2 + 1):
        r = derive('fd', 'face', 'point', 0, list(range(1 - k, k + 2)))
        left = derive('fd', 'face', 'point', 0, list(range(-k, k + 1)))
        for o in range(1 - k, k + 2):
            if left.weights[1 - o] != r.weights[o]:
                raise AssertionError(f'fd biased face k={k}: the left row is not the mirror of the right one')
        for bias, lo, hi in (('br', 1 - k, k + 1), ('bl', -k, k)):
            fw = dict(r.weights) if bias == 'br' else {1 - o: w for o, w in r.weights.items()}
            add(Kernel(f'gideon_fd_face_{bias}{r.order}', 'fd_face_b', 'fd', 'face', bias, 0, k, r.order, lo, hi, 0,
                       dict(r.weights), sum(abs(v) for v in r.weights.values()), fw))
    return ks


# ---------------------------------------------------------------------------------------------------------------------
# cross-check against ADAM's committed tables (parsed, never adopted)
# ---------------------------------------------------------------------------------------------------------------------
NUM = r'[+-]?\d+(?:\.\d*)?_R8P'


def parse_adam_tables(path: Path) -> dict[str, list[Fr]]:
    """Parse ADAM's `real(R8P), parameter :: NAME_Ss(...)=[...]/den` rows into exact fractions."""
    text = path.read_text()
    stmts: list[str] = []
    buf = ''
    for line in text.splitlines():
        code = line.split('!')[0].rstrip()
        if buf:
            buf += ' ' + code.strip().lstrip('&')
        elif 'parameter ::' in code and re.search(r'_S\d\(', code):
            buf = code
        else:
            continue
        if buf.endswith('&'):
            buf = buf[:-1]
            continue
        stmts.append(buf)
        buf = ''
    out: dict[str, list[Fr]] = {}
    for st in stmts:
        m = re.search(r'::\s*(\w+_S\d)\([^)]*\)\s*=\s*\[(.*)\]\s*(?:/\s*(' + NUM + r'))?\s*$', st)
        if not m:
            continue
        name, body, den = m.group(1), m.group(2), m.group(3)
        vals = [Fr(v.replace('_R8P', '').rstrip('.')) for v in re.findall(NUM, body)]
        d = Fr(den.replace('_R8P', '').rstrip('.')) if den else Fr(1)
        out[name] = [v / d for v in vals]
    return out


def cross_check_adam(ks: list[Kernel]) -> list[str]:
    """Compare every ADAM row with the derived one; return the findings (mismatches and coverage)."""
    t = parse_adam_tables(ADAM_REF)
    if not t:
        raise AssertionError(f'no table parsed from {ADAM_REF}')
    by = {(k.family, k.n, k.s, k.bias): k for k in ks}
    found: list[str] = []
    checked = 0

    def cmp(adam: str, mine: list[Fr], vals: list[Fr]) -> None:
        nonlocal checked
        checked += 1
        if mine != vals[: len(mine)] or any(v != 0 for v in vals[len(mine):]):
            found.append(f'MISMATCH {adam}: ADAM {[str(v) for v in vals]} derived {[str(v) for v in mine]}')

    for name, vals in sorted(t.items()):
        m = re.fullmatch(r'(FD|FV|FVR|FVUR|FVUL|DECONV_)(\d|A2P|P2A)?_?(CC|UR|UL)?_S(\d)', name)
        if not m:
            continue
        s = int(m.group(4))
        if all(v == 0 for v in vals):
            continue  # zero rows: absent kernels in GIDEON (audit item C)
        if name.startswith('FD') and name[2].isdigit():
            n = int(name[2])
            k = by.get(('fd_face_c', 0, s, 'c')) if n == 0 else by.get(('fd_node_c', n, s, 'c'))
            if k is None:
                found.append(f'NOT GENERATED {name}')
                continue
            mine = [k.row[i] for i in sorted(k.row)]
            cmp(name, mine, vals)
        elif name.startswith('FV1_CC'):
            k = by[('fv_face_c', 0, s, 'c')]
            cmp(name, [k.row[i] for i in sorted(k.row)], vals)
        elif name.startswith('FVR'):
            n = int(name[3])
            k = by.get(('fv_facedn_c', n, s, 'c'))
            if k is None:
                found.append(f'NOT GENERATED {name}')
                continue
            cmp(name, [k.row[i] for i in sorted(k.row)], vals)
        elif name.startswith('DECONV_A2P') or name.startswith('DECONV_P2A'):
            fam_data = 'fv' if 'A2P' in name else 'fd'
            k = next(x for x in ks if x.family == 'deconv' and x.data == fam_data and x.s == s)
            cmp(name, [k.row[i] for i in sorted(k.row)], vals)
        elif name.startswith('FV1_UR'):  # cells 1..s, face at the left of cell 1: GIDEON's ur face row
            k = by[('fv_face_u', 0, s, 'ur')]
            cmp(name, [k.row[o] for o in range(1, s + 1)], vals)
        elif name.startswith('FV1_UL'):  # stored ascending on cells 1-s..0: GIDEON's ul face row (the mirror)
            k = by[('fv_face_u', 0, s, 'ul')]
            cmp(name, [k.row[1 - c] for c in range(1 - s, 1)], vals)
        elif name.startswith('FVUR'):
            # FVURn(k) multiplies cell 1 - k, the nearest first (left-side stencil): GIDEON's ul face_dn. ADAM's
            # docstring says q(k) = cell c-(k-1) inside the kernel; its callers pass the window ascending, so the
            # coefficient k meets cell c-(k-1) through the reversed access q(s-k+1): the table is right, the
            # docstring describes the window the wrong way round.
            n = int(name[4])
            k = by.get(('fv_facedn_u', n, s, 'ul'))
            if k is None:
                found.append(f'NOT GENERATED {name}')
                continue
            sign = -1 if n % 2 else 1
            cmp(name, [sign * k.row[1 - c] for c in range(0, -s, -1)], vals)
        elif name.startswith('FVUL'):  # FVULn(k) multiplies cell s - k + 1 (right-side stencil): ur face_dn, descending
            n = int(name[4])
            k = by.get(('fv_facedn_u', n, s, 'ur'))
            if k is None:
                found.append(f'NOT GENERATED {name}')
                continue
            cmp(name, [k.row[s - i + 1] for i in range(1, s + 1)], vals)
    found.append(f'checked {checked} ADAM rows')
    return found



# ---------------------------------------------------------------------------------------------------------------------
# orders on a generic smooth function, in 100-digit decimal arithmetic (G3 of the rows, free of round-off)
# ---------------------------------------------------------------------------------------------------------------------
def _sin_cos(x: Decimal) -> tuple[Decimal, Decimal]:
    """sin and cos of a small argument by their Taylor series, at the context precision."""
    s, c, term, k = Decimal(0), Decimal(0), Decimal(1), 0
    eps = Decimal(10) ** (-getcontext().prec - 5)
    while True:
        if k % 4 == 0:
            c += term
        elif k % 4 == 1:
            s += term
        elif k % 4 == 2:
            c -= term
        else:
            s -= term
        k += 1
        term = term * x / k
        if abs(term) < eps and k > 4:
            return s, c


GA, GP, GB, GC, GQ, GX0 = (Decimal('1.3'), Decimal('0.4'), Decimal('0.5'), Decimal('0.7'), Decimal('-0.2'),
                          Decimal('0.37'))


def g_deriv(k: int, x: Decimal) -> Decimal:
    """k-th derivative (k >= -1, -1 the primitive) of g = sin(1.3 x + 0.4) + 0.5 cos(0.7 x - 0.2)."""
    s1, c1 = _sin_cos(GA * x + GP)
    s2, c2 = _sin_cos(GC * x + GQ)
    ds = [s1, c1, -s1, -c1][k % 4]
    dc = [c2, -s2, -c2, s2][k % 4]
    return GA**k * ds + GB * GC**k * dc


def kernel_error(k: Kernel, h: Decimal) -> Decimal:
    """|kernel - target| on g at step h around GX0, every weight exact."""
    geom, tgt = GEOM.get(k.family, ('fvcell', 'point') if k.data == 'fv' else ('fdnode', 'cellavg'))
    half = Decimal('0.5')
    acc = Decimal(0)
    for o, w in k.full.items():
        if geom == 'fdnode':
            d = g_deriv(0, GX0 + o * h)
        elif geom == 'fdface':
            d = g_deriv(0, GX0 + (o - half) * h)
        elif geom == 'fvface':
            d = (g_deriv(-1, GX0 + o * h) - g_deriv(-1, GX0 + (o - 1) * h)) / h
        else:
            d = (g_deriv(-1, GX0 + (o + half) * h) - g_deriv(-1, GX0 + (o - half) * h)) / h
        acc += Decimal(w.numerator) / Decimal(w.denominator) * d
    acc /= h**k.n
    if tgt == 'point':
        t = g_deriv(k.n, GX0)
    else:
        t = (g_deriv(k.n - 1, GX0 + half * h) - g_deriv(k.n - 1, GX0 - half * h)) / h
    return abs(acc - t)


def check_orders(ks: list[Kernel], observed: dict[str, float] | None = None) -> list[str]:
    """Observed order of every row between h = 2^-11 and 2^-12; frozen bound: >= design - 0.25. The observed orders
    are also stored into `observed` (by kernel name) when it is given."""
    getcontext().prec = 100
    out = []
    worst = None
    for k in ks:
        e1, e2 = kernel_error(k, Decimal(2) ** -11), kernel_error(k, Decimal(2) ** -12)
        p = float((e1 / e2).ln() / Decimal(2).ln())
        if observed is not None:
            observed[k.name] = p
        if worst is None or p - k.order < worst[0]:
            worst = (p - k.order, k.name, p)
        if p < k.order - 0.25:
            out.append(f'ORDER {k.name}: observed {p:.3f}, design {k.order}')
    out.append(f'orders: {len(ks)} rows, smallest margin {worst[0]:+.3f} ({worst[1]}, observed {worst[2]:.3f})')
    return out

# ---------------------------------------------------------------------------------------------------------------------
# emission: Fortran sources (generated, never hand-edited)
# ---------------------------------------------------------------------------------------------------------------------
HEADER = '!< Generated by scripts/gen_gideon.py: do not edit; change the generator and run it again.\n'
GEOM = {  # family -> (geometry of the data, target kind) for the test registry
    'fd_face_c': ('fdface', 'point'),
    'fd_node_c': ('fdnode', 'point'),
    'fd_facedn_c': ('fdface', 'point'),
    'fd_face_u': ('fdface', 'point'),
    'fd_face_b': ('fdface', 'point'),
    'fv_face_c': ('fvface', 'point'),
    'fv_facedn_c': ('fvface', 'point'),
    'fv_face_u': ('fvface', 'point'),
    'fv_facedn_u': ('fvface', 'point'),
    'fv_cell_c': ('fvcell', 'cellavg'),
    'fv_cell_u': ('fvcell', 'cellavg'),
}
GEOM_CODE = {'fdnode': 1, 'fdface': 2, 'fvface': 3, 'fvcell': 4}
TARGET_CODE = {'point': 1, 'cellavg': 2}


def lit(v: Fr) -> str:
    """Fortran literal of an exact rational, `p._R_P/q._R_P` (correctly rounded by the compiler's constant folding)."""
    if v.denominator == 1:
        return f'{v.numerator}._R_P'
    return f'{v.numerator}._R_P/{v.denominator}._R_P'


def row_name(k: Kernel) -> str:
    """Name of the stored row a kernel reads."""
    if k.family == 'fd_face_c':
        return f'FD_FACE_S{k.s}'
    if k.family == 'fd_node_c':
        return f'FD_D{k.n}_S{k.s}'
    if k.family == 'fd_facedn_c':
        return f'FD_FACE_D{k.n}_S{k.s}'
    if k.family == 'fd_face_u':
        return f'FD_UFACE_S{k.s}'
    if k.family == 'fd_face_b':
        return f'FD_BFACE_S{k.s}'
    if k.family in ('fv_face_c', 'fv_facedn_c'):
        return f'FV_FACE_S{k.s}' if k.n == 0 else f'FV_FACE_D{k.n}_S{k.s}'
    if k.family == 'fv_cell_c':
        return f'FV_FACE_S{k.s}' if k.n == 1 else f'FV_FACE_D{k.n - 1}_S{k.s}'
    if k.family in ('fv_face_u', 'fv_facedn_u'):
        return f'FV_UFACE_S{k.s}' if k.n == 0 else f'FV_UFACE_D{k.n}_S{k.s}'
    if k.family == 'fv_cell_u':
        return f'FV_UFACE_S{k.s}' if k.n == 1 else f'FV_UFACE_D{k.n - 1}_S{k.s}'
    if k.family == 'deconv':
        return f'FV_POINT_S{k.s}' if k.data == 'fv' else f'FD_AVERAGE_S{k.s}'
    raise ValueError(k.family)


def array_constructor(vals: list[Fr], indent: int) -> str:
    """`[a, b, ...]` wrapped at 120 columns with continuation lines."""
    items = [lit(v) for v in vals]
    lines, cur = [], '['
    for i, it in enumerate(items):
        piece = it + (', ' if i < len(items) - 1 else ']')
        if len(cur) + len(piece) + indent > 116:
            lines.append(cur.rstrip() + ' &')
            cur = ' ' * 4 + piece
        else:
            cur += piece
    lines.append(cur)
    return ('\n' + ' ' * indent).join(lines)


def emit_tables(ks: list[Kernel]) -> str:
    """The rows, one named constant each, stored once (left-upwind and cell kernels read the rows they share)."""
    rows: dict[str, tuple[int, int, list[Fr], str]] = {}
    for k in ks:
        if k.family in ('fv_cell_c', 'fv_cell_u', ) or k.bias in ('ul', 'bl'):
            continue
        name = row_name(k)
        lo = min(k.row)
        vals = [k.row[i] for i in range(lo, max(k.row) + 1)]
        layout = {'c': 'pairs', 'ur': 'cells 1..S, the face nearest first' if k.data == 'fv'
                  else 'points 1..S, the face nearest first', 'br': f'points {lo}..{max(k.row)}, by offset'}[k.bias]
        what = f'{k.name} ({layout})'
        rows[name] = (lo, max(k.row), vals, what)
    out = [HEADER + '!< GIDEON, the stencil rows: exact rationals, one named constant per row, folded by the compiler.',
           'module gideon_tables',
           '!< GIDEON, the stencil rows: exact rationals, one named constant per row, folded by the compiler.',
           '!<',
           '!< Centred rows hold half the stencil in pair form (index m multiplies the pair of offsets m and -m',
           '!< for a node or cell, m and 1-m for a face; index 0 the centre); upwind rows hold the right-side',
           '!< stencil, cells (or points) 1..S, the face nearest first; biased rows hold the right-biased stencil by',
           '!< offset; the left-side kernels read the same rows mirrored.',
           '!< See docs/guide/conventions.md.',
           'use penf, only : R_P', 'implicit none', 'public']
    for name, (lo, hi, vals, what) in rows.items():
        decl = f'real(R_P), parameter :: {name}({lo}:{hi}) = '
        out.append(f'!< {what}')
        out.append(decl + array_constructor(vals, len(decl)))
    out.append('endmodule gideon_tables')
    return '\n'.join(out) + '\n'


def scale(n: int, style: str) -> str:
    """The division by the step: ADAM's FD form (repeated divisions) or FV form (one product), D4."""
    if style == 'fd':
        return '/ds' * n
    return '/ds' if n == 1 else '/(' + '*'.join(['ds'] * n) + ')'


def ordinal(n: int) -> str:
    return {0: 'value', 1: 'first derivative', 2: 'second derivative', 3: 'third derivative',
            4: 'fourth derivative', 5: 'fifth derivative', 6: 'sixth derivative'}[n]


def emit_kernel(k: Kernel) -> str:
    """One kernel: explicit-shape window, the row as a compile-time constant, ADAM's arithmetic shape."""
    rn = row_name(k)
    s = k.s
    has_ds = k.n >= 1
    res = {'face': 'f', 'point': 'v', 'average': 'v'}.get(k.op, 'd')
    args = f'(q, ds, {res})' if has_ds else f'(q, {res})'
    data_doc = 'point values' if k.data == 'fd' else 'cell averages'
    where = {'fd_face_c': 'at the face between nodes 0 and 1', 'fd_node_c': 'at node 0',
             'fd_facedn_c': 'at the face between nodes 0 and 1', 'fd_face_u': 'at the face between nodes 0 and 1',
             'fd_face_b': 'at the face between nodes 0 and 1',
             'fv_face_c': 'at the face between cells 0 and 1', 'fv_facedn_c': 'at the face between cells 0 and 1',
             'fv_face_u': 'at the face between cells 0 and 1', 'fv_facedn_u': 'at the face between cells 0 and 1',
             'fv_cell_c': 'averaged over cell 0', 'fv_cell_u': 'averaged over cell 0'}.get(k.family, '')
    bias_doc = {'c': 'centred', 'ur': 'upwind, right-side stencil', 'ul': 'upwind, left-side stencil',
                'br': 'upwind-biased, one point more on the right', 'bl': 'upwind-biased, one point more on the left'}[k.bias]
    if k.family == 'fd_face_u':
        bias_doc = 'one-sided (extrapolated), ' + ('right-side stencil' if k.bias == 'ur' else 'left-side stencil')
    if k.family == 'deconv':
        title = ('Point value at the centre of cell 0 from cell averages' if k.data == 'fv'
                 else 'Average over cell 0 from point values')
    else:
        title = f'{ordinal(k.n).capitalize()} {where}, {bias_doc}'
    doc = f'   !< {title}, order {k.order}, from {data_doc} q({k.lo}:{k.hi}).'
    if len(doc) > 118:
        doc = f'   !< {title}, order {k.order},\n   !< from {data_doc} q({k.lo}:{k.hi}).'
    L = [f'   pure subroutine {k.name}{args}', doc,
         f'   real(R_P),    intent(in)  :: q({k.lo}:{k.hi}) !< Data by offset (docs/guide/conventions.md).']
    if has_ds:
        L.append('   real(R_P),    intent(in)  :: ds       !< Grid step.')
    L.append(f'   real(R_P),    intent(out) :: {res}        !< Result.')
    fam = k.family
    if fam in ('fv_cell_c', 'fv_cell_u'):
        L.append('   real(R_P)                 :: gl, gr   !< The face rows at the faces -1/2 and +1/2.')
    L.append('   integer(I4P)              :: m        !< Counter.')
    L += ['   !$acc routine seq', '   !$omp declare target']
    pm = '+' if k.parity == 1 else '-'
    if fam in ('fd_face_c', 'fd_facedn_c', 'fv_face_c', 'fv_facedn_c'):
        L += [f'   {res} = 0._R_P', f'   do m=1, {s}', f'      {res} = {res} + {rn}(m)*(q(m) {pm} q(1-m))', '   enddo']
        if has_ds:
            L.append(f'   {res} = {res}{scale(k.n, k.data)}')  # D4: fd repeated divisions, fv one product
    elif fam == 'fd_node_c' or fam == 'deconv':
        if k.parity == 1:
            L.append(f'   {res} = {rn}(0)*q(0)')
        else:
            L.append(f'   {res} = 0._R_P')
        L += [f'   do m=1, {s}', f'      {res} = {res} + {rn}(m)*(q(m) {pm} q(-m))', '   enddo']
        if has_ds:
            L.append(f'   {res} = {res}{scale(k.n, "fd")}')
    elif fam == 'fv_cell_c':
        L += ['   gl = 0._R_P', f'   do m=1, {s}', f'      gl = gl + {rn}(m)*(q(m-1) {pm} q(-m))', '   enddo',
              '   gr = 0._R_P', f'   do m=1, {s}', f'      gr = gr + {rn}(m)*(q(m) {pm} q(1-m))', '   enddo',
              f'   d = (gr-gl){scale(k.n, "fv")}']
    elif fam in ('fv_face_u', 'fv_facedn_u', 'fd_face_u'):
        acc = 'q(m)' if k.bias == 'ur' else 'q(1-m)'
        L += [f'   {res} = 0._R_P', f'   do m=1, {s}', f'      {res} = {res} + {rn}(m)*{acc}', '   enddo']
        if has_ds:
            neg = '-' if (k.bias == 'ul' and k.n % 2) else ''
            L.append(f'   {res} = {neg}{res}{scale(k.n, "fv")}')
    elif fam == 'fd_face_b':
        acc = 'q(m)' if k.bias == 'br' else 'q(1-m)'
        L += [f'   {res} = 0._R_P', f'   do m={1 - s}, {s + 1}', f'      {res} = {res} + {rn}(m)*{acc}', '   enddo']
    elif fam == 'fv_cell_u':
        if k.bias == 'ur':
            L += ['   gl = 0._R_P', f'   do m=1, {s}', f'      gl = gl + {rn}(m)*q(m-1)', '   enddo',
                  '   gr = 0._R_P', f'   do m=1, {s}', f'      gr = gr + {rn}(m)*q(m)', '   enddo',
                  f'   d = (gr-gl){scale(k.n, "fv")}']
        else:  # mirrored: the face +1/2 reads cells 0, -1, ...; the face -1/2 cells -1, -2, ...; rows times (-1)^(n-1)
            first, second = ('gr', 'gl') if (k.n - 1) % 2 == 0 else ('gl', 'gr')
            L += ['   gl = 0._R_P', f'   do m=1, {s}', f'      gl = gl + {rn}(m)*q(-m)', '   enddo',
                  '   gr = 0._R_P', f'   do m=1, {s}', f'      gr = gr + {rn}(m)*q(1-m)', '   enddo',
                  f'   d = ({first}-{second}){scale(k.n, "fv")}']
    else:
        raise ValueError(fam)
    L.append(f'   endsubroutine {k.name}')
    return '\n'.join(L)


MODULES = {  # module -> (families, title)
    'gideon_fd_centred': (('fd_face_c', 'fd_node_c', 'fd_facedn_c'), 'finite differences on point values, centred'),
    'gideon_fd_upwind': (('fd_face_u', 'fd_face_b'), 'finite differences on point values, one-sided and upwind-biased'),
    'gideon_fv_centred': (('fv_face_c', 'fv_facedn_c', 'fv_cell_c'), 'finite volumes on cell averages, centred'),
    'gideon_fv_upwind': (('fv_face_u', 'fv_facedn_u', 'fv_cell_u'), 'finite volumes on cell averages, upwind'),
    'gideon_deconvolution': (('deconv',), 'cell averages and point values, one from the other'),
}


def emit_module(mod: str, ks: list[Kernel]) -> str:
    fams, title = MODULES[mod]
    mine = [k for k in ks if k.family in fams]
    out = [HEADER + f'!< GIDEON, {title}.', f'module {mod}', f'!< GIDEON, {title}.', '!<',
           '!< One pure, device-callable kernel per operator and order: the window an explicit-shape dummy by offset,',
           '!< the row a compile-time constant of gideon_tables, ADAM\'s arithmetic shape (pair form, division by the',
           '!< step inside). Conventions and names: docs/guide/conventions.md.',
           'use penf, only : I4P, R_P', 'use gideon_tables', 'implicit none', 'private']
    out += [f'public :: {k.name}' for k in mine]
    out.append('contains')
    out += [emit_kernel(k) + '\n' for k in mine]
    out.append(f'endmodule {mod}')
    return '\n'.join(out) + '\n'


def leading_error(k: Kernel) -> Fr:
    """Exact error of the operator on x^(K+1) at h = 1 (K its exact degree): non-zero by the definition of K."""
    geom, tgt = GEOM.get(k.family, ('fvcell', 'point') if k.data == 'fv' else ('fdnode', 'cellavg'))
    data, place = {'fdnode': ('fd', 'node'), 'fdface': ('fd', 'face'), 'fvface': ('fv', 'face'),
                   'fvcell': ('fv', 'cell')}[geom]
    j = k.order + k.n
    e = sum(w * sampler(data, place, o)(j) for o, w in k.full.items()) - target(tgt, k.n)(j)
    if e == 0:
        raise AssertionError(f'{k.name}: exact beyond its degree')
    return e


def dispatcher_name(k: Kernel) -> str:
    """The run-time order dispatcher of a kernel: its name without the order."""
    return re.sub(r'\d+$', '', k.name)


def _fit(pairs: list[tuple[int, int]]) -> str | None:
    """`order/2 + c` or `order + c` matching every (order, value) pair, or None."""
    for a in (2, 1):
        cs = {v - o // a for o, v in pairs}
        if len(cs) == 1:
            c = cs.pop()
            return ('order/2' if a == 2 else 'order') + (f'+{c}' if c > 0 else f'{c}' if c < 0 else '')
    return None


def _eval(expr: str, order: int) -> int:
    """Evaluate a bound expression as Fortran does (integer division truncates; orders are positive)."""
    return eval(expr.replace('order/2', f'({order}//2)').replace('order', str(order)))  # noqa: S307


def window_expr(group: list[Kernel]) -> tuple[str, str]:
    """Window bounds of a dispatcher as expressions of `order`, verified against every kernel of the group."""
    if len({k.hi for k in group}) == 1:  # left-side upwind: the upper bound fixed, the lower one grows
        hi = str(group[0].hi)
        f = _fit([(k.order, -k.lo) for k in group])
        lo = f'-({f})' if f else None
    else:
        hi = _fit([(k.order, k.hi) for k in group])
        if hi is None:
            lo = None
        elif all(k.lo == -k.hi for k in group):
            lo = f'-({hi})'
        elif all(k.lo == 1 - k.hi for k in group):
            lo = f'1-({hi})'
        elif len({k.lo for k in group}) == 1:
            lo = str(group[0].lo)
        elif len({k.lo + k.hi for k in group}) == 1:  # lo = c - hi (the right-biased windows: c = 2)
            lo = f'{group[0].lo + group[0].hi}-({hi})'
        else:
            lo = None
    if hi is None or lo is None:
        raise AssertionError(f'no window formula for {dispatcher_name(group[0])}')
    for k in group:
        if (_eval(lo, k.order), _eval(hi, k.order)) != (k.lo, k.hi):
            raise AssertionError(f'{k.name}: window {k.lo}:{k.hi} is not {lo}:{hi}')
    return lo, hi


def groups(ks: list[Kernel]) -> dict[str, list[Kernel]]:
    out: dict[str, list[Kernel]] = {}
    for k in ks:
        out.setdefault(dispatcher_name(k), []).append(k)
    return out


def bad_order(group: list[Kernel]) -> int:
    """An order with no kernel whose window stays narrow (the dispatch test's NaN case)."""
    return max(k.order for k in group) + (2 if group[0].bias == 'c' else 1)


def emit_dispatch(ks: list[Kernel]) -> str:
    """Run-time order dispatchers: one select case per operator, the whole window passed to the kernel."""
    title = '!< GIDEON, the run-time order dispatchers: the order chosen when the program runs.'
    out = [HEADER + title, 'module gideon_dispatch', title, '!<',
           '!< `gideon_<data>_<op>_<bias>(order, q, [ds,] r)` calls `gideon_<data>_<op>_<bias><order>`: one select',
           '!< case per call, uniform across cells when the order is (no divergence on a GPU). The window bounds are',
           '!< expressions of order, so the whole window is passed to the kernel. An order without a kernel returns',
           '!< GIDEON_NAN, a quiet NaN (error stop is not device code; a silent 0 is what ADAM returned).',
           'use penf, only : I4P, R_P', 'use gideon_kind_parameters, only : GIDEON_NAN',
           'use gideon_fd_centred', 'use gideon_fd_upwind', 'use gideon_fv_centred', 'use gideon_fv_upwind',
           'use gideon_deconvolution',
           'implicit none', 'private']
    gs = groups(ks)
    out += [f'public :: {name}' for name in gs]
    out.append('contains')
    for name, group in gs.items():
        lo, hi = window_expr(group)
        k0 = group[0]
        has_ds = k0.n >= 1
        res = {'face': 'f', 'point': 'v', 'average': 'v'}.get(k0.op, 'd')
        args = f'(order, q, ds, {res})' if has_ds else f'(order, q, {res})'
        orders = ', '.join(str(k.order) for k in group)
        out += [f'   pure subroutine {name}{args}',
                f'   !< Run-time order dispatch to {name}<order>, order in {{{orders}}}.',
                '   integer(I4P), intent(in)  :: order !< Order of accuracy.',
                f'   real(R_P),    intent(in)  :: q({lo}:{hi}) !< Data by offset, the window of that order.']
        if has_ds:
            out.append('   real(R_P),    intent(in)  :: ds    !< Grid step.')
        out += [f'   real(R_P),    intent(out) :: {res}     !< Result; GIDEON_NAN for an order without a kernel.',
                '   !$acc routine seq', '   !$omp declare target', '   select case(order)']
        for k in group:
            call = f'call {k.name}(q, ds, {res})' if has_ds else f'call {k.name}(q, {res})'
            out.append(f'   case({k.order}) ; {call}')
        out += [f'   case default ; {res} = GIDEON_NAN', '   endselect', f'   endsubroutine {name}', '']
    out.append('endmodule gideon_dispatch')
    return '\n'.join(out) + '\n'


OPS3D = ('gradient', 'divergence', 'curl', 'laplacian')


def kernels_3d(ks: list[Kernel]) -> list[tuple[str, str, str, int, int]]:
    """(name, data, op, order, s) of the 3-D operators: centred, over the 1-D d1 (d2 for the Laplacian) kernels."""
    out = []
    for data, fam in (('fd', 'fd_node_c'), ('fv', 'fv_cell_c')):
        for op in OPS3D:
            n = 2 if op == 'laplacian' else 1
            for k in ks:
                if k.family == fam and k.n == n:
                    out.append((f'gideon_{data}_{op}_c{k.order}', data, op, k.order, k.s))
    return out


def emit_operators_3d(ks: list[Kernel]) -> str:
    """Gradient, divergence, curl, Laplacian over directional stencils (ADAM's FNL signature and sum order)."""
    title = '!< GIDEON, the 3-D operators: gradient, divergence, curl and Laplacian over directional stencils.'
    out = [HEADER + title, 'module gideon_operators_3d', title, '!<',
           '!< The caller gathers one contiguous stencil per direction and component (ADAM\'s FNL signature: no',
           '!< non-contiguous section inside a device routine), cell or node 0 the point of the result: `qx(-S:S)`',
           '!< the field along x, `vy_x(-S:S)` the y component along x. Each operator is the 1-D kernels summed in',
           '!< ADAM\'s order (x, y, z), so it reproduces ADAM\'s value bitwise; the run-time order dispatchers follow.',
           'use penf, only : I4P, R_P', 'use gideon_kind_parameters, only : GIDEON_NAN', 'use gideon_fd_centred',
           'use gideon_fv_centred', 'implicit none', 'private']
    k3 = kernels_3d(ks)
    names = [x[0] for x in k3]
    disp = sorted({re.sub(r'\d+$', '', nm) for nm in names})
    out += [f'public :: {nm}' for nm in names + disp]
    out.append('contains')
    for name, data, op, order, s in k3:
        w = f'-{s}:{s}'
        d1 = f'gideon_{data}_d1_c{order}'
        d2 = f'gideon_{data}_d2_c{order}'
        what = 'point values' if data == 'fd' else 'cell averages'
        if op == 'gradient':
            out += [f'   pure subroutine {name}(qx, qy, qz, dxyz, g)',
                    f'   !< Gradient at point 0, centred, order {order}, from {what} along x, y, z.',
                    f'   real(R_P), intent(in)  :: qx({w}), qy({w}), qz({w}) !< The field along x, y, z.',
                    '   real(R_P), intent(in)  :: dxyz(3) !< Steps.',
                    '   real(R_P), intent(out) :: g(3)    !< Gradient.',
                    '   !$acc routine seq', '   !$omp declare target',
                    f'   call {d1}(qx, dxyz(1), g(1))', f'   call {d1}(qy, dxyz(2), g(2))',
                    f'   call {d1}(qz, dxyz(3), g(3))']
        elif op == 'divergence':
            out += [f'   pure subroutine {name}(vx_x, vy_y, vz_z, dxyz, div)',
                    f'   !< Divergence at point 0, centred, order {order}, from {what}.',
                    f'   real(R_P), intent(in)  :: vx_x({w}), vy_y({w}), vz_z({w}) !< Each component along its axis.',
                    '   real(R_P), intent(in)  :: dxyz(3) !< Steps.',
                    '   real(R_P), intent(out) :: div     !< Divergence.',
                    '   real(R_P)              :: dx, dy, dz',
                    '   !$acc routine seq', '   !$omp declare target',
                    f'   call {d1}(vx_x, dxyz(1), dx)', f'   call {d1}(vy_y, dxyz(2), dy)',
                    f'   call {d1}(vz_z, dxyz(3), dz)', '   div = dx + dy + dz']
        elif op == 'curl':
            out += [f'   pure subroutine {name}(vy_x, vz_x, vx_y, vz_y, vx_z, vy_z, dxyz, c)',
                    f'   !< Curl at point 0, centred, order {order}, from {what}: the six off-axis component stencils.',
                    f'   real(R_P), intent(in)  :: vy_x({w}), vz_x({w}) !< v_y, v_z along x.',
                    f'   real(R_P), intent(in)  :: vx_y({w}), vz_y({w}) !< v_x, v_z along y.',
                    f'   real(R_P), intent(in)  :: vx_z({w}), vy_z({w}) !< v_x, v_y along z.',
                    '   real(R_P), intent(in)  :: dxyz(3) !< Steps.',
                    '   real(R_P), intent(out) :: c(3)    !< Curl.',
                    '   real(R_P)              :: dvx_dy, dvx_dz, dvy_dx, dvy_dz, dvz_dx, dvz_dy',
                    '   !$acc routine seq', '   !$omp declare target',
                    f'   call {d1}(vx_y, dxyz(2), dvx_dy)', f'   call {d1}(vx_z, dxyz(3), dvx_dz)',
                    f'   call {d1}(vy_x, dxyz(1), dvy_dx)', f'   call {d1}(vy_z, dxyz(3), dvy_dz)',
                    f'   call {d1}(vz_x, dxyz(1), dvz_dx)', f'   call {d1}(vz_y, dxyz(2), dvz_dy)',
                    '   c(1) = dvz_dy - dvy_dz', '   c(2) = dvx_dz - dvz_dx', '   c(3) = dvy_dx - dvx_dy']
        else:
            out += [f'   pure subroutine {name}(qx, qy, qz, dxyz, lap)',
                    f'   !< Laplacian at point 0, centred, order {order}, from {what} along x, y, z.',
                    f'   real(R_P), intent(in)  :: qx({w}), qy({w}), qz({w}) !< The field along x, y, z.',
                    '   real(R_P), intent(in)  :: dxyz(3) !< Steps.',
                    '   real(R_P), intent(out) :: lap     !< Laplacian.',
                    '   real(R_P)              :: dx, dy, dz',
                    '   !$acc routine seq', '   !$omp declare target',
                    f'   call {d2}(qx, dxyz(1), dx)', f'   call {d2}(qy, dxyz(2), dy)',
                    f'   call {d2}(qz, dxyz(3), dz)', '   lap = dx + dy + dz']
        out += [f'   endsubroutine {name}', '']
    for dn in disp:
        group = [x for x in k3 if re.sub(r'\d+$', '', x[0]) == dn]
        op = group[0][2]
        hi = 'order/2' if op != 'laplacian' else 'order/2'
        w = f'-({hi}):{hi}'
        for x in group:
            if x[4] != x[3] // 2:
                raise AssertionError(f'{x[0]}: half-stencil is not order/2')
        sig = {'gradient': ('(order, qx, qy, qz, dxyz, g)', ['qx', 'qy', 'qz'], 'g(3)', 'g = GIDEON_NAN'),
               'divergence': ('(order, vx_x, vy_y, vz_z, dxyz, div)', ['vx_x', 'vy_y', 'vz_z'], 'div',
                              'div = GIDEON_NAN'),
               'curl': ('(order, vy_x, vz_x, vx_y, vz_y, vx_z, vy_z, dxyz, c)',
                        ['vy_x', 'vz_x', 'vx_y', 'vz_y', 'vx_z', 'vy_z'], 'c(3)', 'c = GIDEON_NAN'),
               'laplacian': ('(order, qx, qy, qz, dxyz, lap)', ['qx', 'qy', 'qz'], 'lap', 'lap = GIDEON_NAN')}[op]
        args, ins, res, nan = sig
        out += [f'   pure subroutine {dn}{args}',
                f'   !< Run-time order dispatch to {dn}<order>, order in {{{", ".join(str(x[3]) for x in group)}}}.',
                '   integer(I4P), intent(in)  :: order !< Order of accuracy.',
                *[f'   real(R_P),    intent(in)  :: {a}({w}) !< Stencil of that order.' for a in ins],
                '   real(R_P),    intent(in)  :: dxyz(3) !< Steps.',
                f'   real(R_P),    intent(out) :: {res} !< Result; GIDEON_NAN for an order without a kernel.',
                '   !$acc routine seq', '   !$omp declare target', '   select case(order)']
        call_args = args.replace('(order, ', '(')
        for x in group:
            out.append(f'   case({x[3]}) ; call {x[0]}{call_args}')
        out += [f'   case default ; {nan}', '   endselect', f'   endsubroutine {dn}', '']
    out.append('endmodule gideon_operators_3d')
    return '\n'.join(out) + '\n'


def emit_registry(ks: list[Kernel]) -> str:
    """Test registry: metadata of every kernel and one routine that applies kernel `id` to a wide window."""
    nk = len(ks)
    meta = []
    for i, k in enumerate(ks, 1):
        geom, tgt = GEOM.get(k.family, ('fvcell', 'point') if k.data == 'fv' else ('fdnode', 'cellavg'))
        flux = 1 if k.family in ('fv_cell_c', 'fv_cell_u') or (k.family == 'fd_node_c' and k.parity == -1) else 0
        flux = 1 if k.family == 'fd_facedn_c' and k.parity == -1 else flux  # antisymmetric face pairs
        meta.append((i, k, GEOM_CODE[geom], TARGET_CODE[tgt], flux))
    title = '!< GIDEON tests, the registry of the kernels: their metadata and one routine that applies any of them.'
    out = [HEADER + title, 'module gideon_test_registry', title,
           '!<',
           '!< geom: 1 point values at x = o (node 0), 2 point values at x = o - 1/2 (face at 0), 3 cell o on [o-1, o]',
           '!< (face at 0), 4 cell o on [o-1/2, o+1/2] (cell 0). target: 1 the n-th derivative at x = 0, 2 its average',
           '!< over cell 0. zero_const: the kernel returns exactly 0 on constant data (pair or flux-difference form).',
           'use penf, only : I4P, R_P', 'use gideon', 'implicit none', 'private',
           'public :: NK, WMAX, kname, kdata, kbias, kgeom, ktarget, kn, korder, klo, khi, kdegree, ksumabs',
           'public :: kzero_const, kerror',
           'public :: kmirror, apply, apply_dev, apply_dispatch, apply_dispatch_bad',
           f'integer(I4P), parameter :: NK = {nk} !< Number of kernels.',
           'integer(I4P), parameter :: WMAX = 16 !< Half-width of the windows the tests pass.']

    def arr(name: str, typ: str, vals: list[str], ln: str = '') -> None:
        decl = f'{typ}, parameter :: {name}(NK) = '
        lines, cur = [], '['
        for j, v in enumerate(vals):
            piece = v + (', ' if j < len(vals) - 1 else ']')
            if len(cur) + len(piece) + len(decl) > 112:
                lines.append(cur.rstrip() + ' &')
                cur = '   ' + piece
            else:
                cur += piece
        lines.append(cur)
        out.append(decl + ('\n' + ' ' * len(decl)).join(lines) + ln)

    arr('kname', 'character(len=24)', [f"'{k.name.removeprefix('gideon_'):24s}'" for _, k, *_ in meta])
    arr('kdata', 'character(len=2)', [f"'{k.data}'" for _, k, *_ in meta])
    arr('kbias', 'character(len=2)', [f"'{k.bias:2s}'" for _, k, *_ in meta])
    arr('kgeom', 'integer(I4P)', [str(g) for _, _, g, _, _ in meta])
    arr('ktarget', 'integer(I4P)', [str(t) for _, _, _, t, _ in meta])
    arr('kn', 'integer(I4P)', [str(k.n) for _, k, *_ in meta])
    arr('korder', 'integer(I4P)', [str(k.order) for _, k, *_ in meta])
    arr('klo', 'integer(I4P)', [str(k.lo) for _, k, *_ in meta])
    arr('khi', 'integer(I4P)', [str(k.hi) for _, k, *_ in meta])
    arr('kdegree', 'integer(I4P)', [str(k.order + k.n - 1) for _, k, *_ in meta])
    arr('ksumabs', 'real(R_P)', [f'{float(k.sumabs):.17e}'.replace('e', '_R_P*1e') + '_R_P' for _, k, *_ in meta])
    arr('kerror', 'real(R_P)', [lit(leading_error(k)) for _, k, *_ in meta])  # exact: R16P compares at 34 digits
    arr('kzero_const', 'logical', ['.true.' if z else '.false.' for *_, z in meta])
    # mirror partner of each left-side kernel, ul -> ur and bl -> br (0 elsewhere)
    idx = {k.name: i for i, k, *_ in meta}
    partner = {'ul': ('_ul', '_ur'), 'bl': ('_bl', '_br')}
    arr('kmirror', 'integer(I4P)', [str(idx.get(k.name.replace(*partner[k.bias]), 0)) if k.bias in partner else '0'
                                    for _, k, *_ in meta])
    out += ['contains',
            '   subroutine apply(id, q, ds, r)',
            '   !< Apply kernel `id` to its window of q (offsets klo(id)..khi(id)).',
            '   integer(I4P), intent(in)  :: id         !< Kernel.',
            '   real(R_P),    intent(in)  :: q(-WMAX:WMAX) !< Data by offset.',
            '   real(R_P),    intent(in)  :: ds         !< Grid step (ignored by the kernels without one).',
            '   real(R_P),    intent(out) :: r          !< Result.',
            '   select case(id)']
    for i, k, *_ in meta:
        call = f'call {k.name}(q({k.lo}:{k.hi}), ds, r)' if k.n >= 1 else f'call {k.name}(q({k.lo}:{k.hi}), r)'
        out.append(f'   case({i}) ; {call}')
    out += ["   case default ; error stop 'apply: no such kernel'", '   endselect', '   endsubroutine apply', '',
            '   pure subroutine apply_dev(id, q, ds, r)',
            '   !< apply, device-callable: the same calls, GIDEON_NAN for an unknown id (no error stop on a device).',
            '   integer(I4P), intent(in)  :: id', '   real(R_P),    intent(in)  :: q(-WMAX:WMAX)',
            '   real(R_P),    intent(in)  :: ds', '   real(R_P),    intent(out) :: r',
            '   !$acc routine seq', '   !$omp declare target', '   select case(id)']
    for i, k, *_ in meta:
        call = f'call {k.name}(q({k.lo}:{k.hi}), ds, r)' if k.n >= 1 else f'call {k.name}(q({k.lo}:{k.hi}), r)'
        out.append(f'   case({i}) ; {call}')
    out += ['   case default ; r = GIDEON_NAN', '   endselect', '   endsubroutine apply_dev', '']
    gs = groups(ks)
    for bad in (False, True):
        sub = 'apply_dispatch_bad' if bad else 'apply_dispatch'
        what = ('its dispatcher with an order that has no kernel (window of that order)' if bad
                else 'its run-time order dispatcher')
        out += [f'   subroutine {sub}(id, q, ds, r)', f'   !< Apply kernel `id` through {what}.',
                '   integer(I4P), intent(in)  :: id', '   real(R_P),    intent(in)  :: q(-WMAX:WMAX)',
                '   real(R_P),    intent(in)  :: ds', '   real(R_P),    intent(out) :: r', '   select case(id)']
        for i, k, *_ in meta:
            name = dispatcher_name(k)
            lo, hi = window_expr(gs[name])
            order = bad_order(gs[name]) if bad else k.order
            sec = f'q({_eval(lo, order)}:{_eval(hi, order)})'
            call = (f'call {name}({order}, {sec}, ds, r)' if k.n >= 1 else f'call {name}({order}, {sec}, r)')
            out.append(f'   case({i}) ; {call}')
        out += ["   case default ; error stop 'no such kernel'", '   endselect', f'   endsubroutine {sub}', '']
    out.append('endmodule gideon_test_registry')
    return '\n'.join(out) + '\n'


def parse_emitted_tables(text: str) -> dict[str, list[Fr]]:
    """Read the emitted rows back (G1: the text written is the row derived)."""
    out: dict[str, list[Fr]] = {}
    stmts, buf = [], ''
    for line in text.splitlines():
        code = line.split('!')[0].rstrip()
        if buf:
            buf += ' ' + code.strip()
        elif 'parameter ::' in code:
            buf = code
        else:
            continue
        if buf.endswith('&'):
            buf = buf[:-1]
            continue
        stmts.append(buf)
        buf = ''
    for st in stmts:
        m = re.search(r'::\s*(\w+)\(\s*-?\d+:\d+\)\s*=\s*\[(.*)\]', st)
        vals = []
        for item in m.group(2).split(','):
            num = re.fullmatch(r'\s*(-?\d+)\._R_P(?:/(\d+)\._R_P)?\s*', item)
            vals.append(Fr(int(num.group(1)), int(num.group(2) or 1)))
        out[m.group(1)] = vals
    return out


def emit_docs_table(ks: list[Kernel]) -> str:
    """The kernel table of the operators page (generated, never typed: DOC5)."""
    out = ['<!-- Generated by scripts/gen_gideon.py: do not edit. -->', '',
           '| Kernel | Data | Result | Bias | Order | Window | Run-time dispatcher |',
           '|---|---|---|---|---|---|---|']
    what = {'face': 'value at the face', 'point': 'point value at the cell centre', 'average': 'cell average'}
    for k in ks:
        res = what.get(k.op) or (f'd{k.n} at the face' if k.op.startswith('face_') else
                                 (f'd{k.n} at the node' if k.data == 'fd' else f'd{k.n} averaged over the cell'))
        data = 'point values' if k.data == 'fd' else 'cell averages'
        bias = {'c': 'centred', 'ur': 'right side', 'ul': 'left side', 'br': 'biased right', 'bl': 'biased left'}[k.bias]
        out.append(f'| `{k.name}` | {data} | {res} | {bias} | {k.order} | `q({k.lo}:{k.hi})` | '
                   f'`{dispatcher_name(k)}` |')
    for name, data, op, order, s in kernels_3d(ks):
        what3 = 'point values' if data == 'fd' else 'cell averages'
        out.append(f'| `{name}` | {what3}, one stencil per axis | {op} | centred | {order} | `({-s}:{s})` per axis | '
                   f'`{re.sub(r"[0-9]+$", "", name)}` |')
    return '\n'.join(out) + '\n'


def write_sources(ks: list[Kernel]) -> list[str]:
    files = {ROOT / 'src/lib/gideon_tables.F90': emit_tables(ks)}
    for mod in MODULES:
        files[ROOT / f'src/lib/{mod}.F90'] = emit_module(mod, ks)
    files[ROOT / 'src/lib/gideon_dispatch.F90'] = emit_dispatch(ks)
    files[ROOT / 'src/lib/gideon_operators_3d.F90'] = emit_operators_3d(ks)
    files[ROOT / 'src/tests/gideon_test_registry.F90'] = emit_registry(ks)
    files[ROOT / 'docs/guide/_operators-table.md'] = emit_docs_table(ks)
    back = parse_emitted_tables(files[ROOT / 'src/lib/gideon_tables.F90'])
    for k in ks:
        if k.family in ('fv_cell_c', 'fv_cell_u') or k.bias in ('ul', 'bl'):
            continue
        lo = min(k.row)
        if back[row_name(k)] != [k.row[i] for i in range(lo, max(k.row) + 1)]:
            raise AssertionError(f'emitted row {row_name(k)} is not the derived one')
    written = []
    for path, text in files.items():
        if not path.exists() or path.read_text() != text:
            path.write_text(text)
            written.append(str(path.relative_to(ROOT)))
    return written


def emit_json(ks: list[Kernel], findings: list[str], observed: dict[str, float]) -> str:
    """The derived rows as data for the documentation (docs/scripts): every kernel with its window, its weights by
    offset as exact fractions (the full operator: left-side rows mirrored, cell derivatives as flux differences),
    its design order, its order observed at 100 digits, its exact leading error; the 3-D operators; the findings."""
    out = {
        'generator': 'scripts/gen_gideon.py',
        'findings': findings,
        'kernels': [{
            'name': k.name, 'family': k.family, 'data': k.data, 'op': k.op, 'bias': k.bias, 'n': k.n, 's': k.s,
            'order': k.order, 'lo': k.lo, 'hi': k.hi, 'dispatcher': dispatcher_name(k),
            'weights': {str(o): str(w) for o, w in sorted(k.full.items())}, 'sumabs': str(k.sumabs),
            'leading_error': str(leading_error(k)), 'observed_order': observed[k.name],
        } for k in ks],
        'operators_3d': [{'name': name, 'data': data, 'op': op, 'order': order, 's': s}
                         for name, data, op, order, s in kernels_3d(ks)],
        'dispatchers': [{'name': name, 'orders': [k.order for k in g], 'window': list(window_expr(g)),
                         'nan_order': bad_order(g)} for name, g in groups(ks).items()],
    }
    return json.dumps(out, indent=1) + '\n'


# ---------------------------------------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--check', action='store_true', help='derive and verify only; write nothing')
    ap.add_argument('--list', action='store_true', help='list the kernels')
    ap.add_argument('--json', type=Path, metavar='FILE', help='write the rows, windows, orders and findings as JSON')
    args = ap.parse_args()
    ks = build()
    observed: dict[str, float] = {}
    findings = cross_check_adam(ks) + check_orders(ks, observed)
    if args.list:
        for k in ks:
            print(f'{k.name:28s} {k.family:12s} n={k.n} s={k.s} order={k.order} window={k.lo}:{k.hi}')
    for f in findings:
        print(f)
    print(f'{len(ks)} kernels derived')
    bad = [f for f in findings if f.startswith(('MISMATCH', 'NOT GENERATED', 'ORDER'))]
    if bad:
        return 1
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(emit_json(ks, findings, observed))
        print(f'written {args.json}')
    if not args.check:
        for w in write_sources(ks):
            print(f'written {w}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
