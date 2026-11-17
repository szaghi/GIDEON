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
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
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
    family: str  # fd_face_c, fd_node_c, fv_face_c, fv_facedn_c, fv_cell_c, fv_face_u, fv_facedn_u, fv_cell_u, deconv
    data: str
    op: str  # face, dN, face_dN, point, average
    bias: str  # c, ur, ul
    n: int
    s: int  # half-stencil (centred) or width (upwind)
    order: int
    lo: int
    hi: int
    parity: int = 0  # centred: +1 symmetric pairs, -1 antisymmetric pairs
    row: dict[int, Fr] = field(default_factory=dict)  # stored row (centred: pair half-row incl. c0; upwind: ur row)
    sumabs: Fr = Fr(0)  # sum of |weights| of the full operator (round-off scale of the tests)


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
                   sum(abs(v) for v in r.weights.values())))
    for n in range(1, NDER_MAX + 1):
        for s in S_CENTRED:
            if 2 * s + 1 <= n:
                continue
            r = derive('fd', 'node', 'point', n, list(range(-s, s + 1)))
            if r.is_zero() or r.order < 1:
                continue
            par, half = check_pairs(r, 'node')
            add(Kernel(f'gideon_fd_d{n}_c{r.order}', 'fd_node_c', 'fd', f'd{n}', 'c', n, s, r.order, -s, s, par, half,
                       sum(abs(v) for v in r.weights.values())))

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
                       sum(abs(v) for v in r.weights.values())))
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
                       sum(abs(v) for v in full.values())))

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
                add(Kernel(f'gideon_fv_{op}_{bias}{r.order}', fam, 'fv', op, bias, n, s, r.order, lo, hi, 0,
                           dict(r.weights), sum(abs(v) for v in r.weights.values())))
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
                add(Kernel(f'gideon_fv_d{n}_{bias}{order}', 'fv_cell_u', 'fv', f'd{n}', bias, n, s, order, lo, hi, 0,
                           dict(g.weights), sum(abs(v) for v in full.values())))

    # deconvolution: cell average -> point value at the centre, point values -> cell average
    for s in S_CENTRED:
        r = derive('fv', 'cell', 'point', 0, list(range(-s, s + 1)))
        par, half = check_pairs(r, 'node')
        add(Kernel(f'gideon_fv_point_c{r.order}', 'deconv', 'fv', 'point', 'c', 0, s, r.order, -s, s, par, half,
                   sum(abs(v) for v in r.weights.values())))
        r = derive('fd', 'cell', 'cellavg', 0, list(range(-s, s + 1)))
        par, half = check_pairs(r, 'node')
        add(Kernel(f'gideon_fd_average_c{r.order}', 'deconv', 'fd', 'average', 'c', 0, s, r.order, -s, s, par, half,
                   sum(abs(v) for v in r.weights.values())))
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
# main
# ---------------------------------------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--check', action='store_true', help='derive and verify only; write nothing')
    ap.add_argument('--list', action='store_true', help='list the kernels')
    args = ap.parse_args()
    ks = build()
    findings = cross_check_adam(ks)
    if args.list:
        for k in ks:
            print(f'{k.name:28s} {k.family:12s} n={k.n} s={k.s} order={k.order} window={k.lo}:{k.hi}')
    for f in findings:
        print(f)
    print(f'{len(ks)} kernels derived')
    bad = [f for f in findings if f.startswith(('MISMATCH', 'NOT GENERATED'))]
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
