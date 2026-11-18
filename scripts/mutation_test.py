#!/usr/bin/env python3
"""Mutation testing of GIDEON's test suite: does every test-worthy defect make some test fail?

Each mutant changes one thing in a scratch copy of the tree (a coefficient of a generated row, the sign of a pair,
a window offset, the division by the step, the mirroring of a left-side kernel, the order of a flux difference, a
dispatcher's case or its NaN, the wiring of a 3-D operator; #2: a dissipation constant or sign, a filter row, a compact
band or right side, a filter's alpha_f term, a solver's pivot, read, border or mask), the tests are rebuilt
(--force-compile) and run with
scripts/run_tests.sh, and the mutant is *killed* when at least one test fails. A surviving mutant is a defect the
suite cannot see: the report lists them. SETH's harness (scripts/mutation_test.py there), GIDEON's mutants.

Usage (from the repository root, with fobis and the compilers of the mode on PATH, PENF fetched):

    python3 scripts/mutation_test.py [--mode tests-gnu] [--workdir DIR] [--only NAME ...]

The working tree is never modified: the copy goes to --workdir (default: a new temporary directory). Exit status 0
when every mutant is killed, 1 otherwise. The mutants edit the generated Fortran, not the generator: the tests are
what is measured, and they never run the generator.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# (name, file, anchor, old, new): replace the first `old` after `anchor` (the anchor itself when old == anchor). Anchor
# in code, never in a comment: a mutated comment is an equivalent mutant.
T = 'src/lib/gideon_tables.F90'
FDC = 'src/lib/gideon_fd_centred.F90'
FDU = 'src/lib/gideon_fd_upwind.F90'
FVC = 'src/lib/gideon_fv_centred.F90'
FVU = 'src/lib/gideon_fv_upwind.F90'
DEC = 'src/lib/gideon_deconvolution.F90'
DIS = 'src/lib/gideon_dispatch.F90'
O3D = 'src/lib/gideon_operators_3d.F90'
DSS = 'src/lib/gideon_dissipation.F90'
FIL = 'src/lib/gideon_filters.F90'
LS = 'src/lib/gideon_line_solvers.F90'
CMP = 'src/lib/gideon_compact.F90'
MUTANTS = [
    # rows: a coefficient of each family, one subtle (the last, smallest one of the widest FD d1)
    ('row-fd-d1-c6', T, 'FD_D1_S3(1:3)', '3._R_P/4._R_P', '0.7501_R_P'),
    ('row-fd-d1-c10-tail', T, 'FD_D1_S5(1:5)', '1._R_P/1260._R_P', '1._R_P/1261._R_P'),
    ('row-fv-face-c4', T, 'FV_FACE_S2(1:2)', '7._R_P/12._R_P', '0.5834_R_P'),
    ('row-fv-uface-3', T, 'FV_UFACE_S3(1:3)', '-7._R_P/6._R_P', '-1.1667_R_P'),
    ('row-fv-point-c6', T, 'FV_POINT_S2(0:2)', '3._R_P/640._R_P', '3._R_P/641._R_P'),
    # pair signs and windows of the centred kernels
    ('pair-sign-fd-d2-c4', FDC, 'pure subroutine gideon_fd_d2_c4', '(q(m) + q(-m))', '(q(m) - q(-m))'),
    ('pair-offset-fd-face-c8', FDC, 'pure subroutine gideon_fd_face_c8', 'q(m) + q(1-m)', 'q(m) + q(-m)'),
    ('pair-sign-fd-average-c6', DEC, 'pure subroutine gideon_fd_average_c6', '(q(m) + q(-m))', '(q(m) - q(-m))'),
    ('flux-face-offset-fv-d1-c6', FVC, 'pure subroutine gideon_fv_d1_c6', '(q(m-1) + q(-m))', '(q(m) + q(-m))'),
    # the division by the step
    ('scale-fd-d3-c4', FDC, 'pure subroutine gideon_fd_d3_c4', 'd = d/ds/ds/ds', 'd = d/ds/ds'),
    ('scale-fv-d2-c4', FVC, 'pure subroutine gideon_fv_d2_c4', '/(ds*ds)', '/ds'),
    ('scale-fv-d6-ur2', FVU, 'pure subroutine gideon_fv_d6_ur2', '/(ds*ds*ds*ds*ds*ds)', '/(ds*ds*ds*ds*ds)'),
    # upwind: the mirror, its sign, the flux difference
    ('mirror-fv-face-ul3', FVU, 'pure subroutine gideon_fv_face_ul3', 'q(1-m)', 'q(-m)'),
    ('mirror-sign-fv-face-d1-ul2', FVU, 'pure subroutine gideon_fv_face_d1_ul2', 'd = -d/ds', 'd = d/ds'),
    ('difference-order-fv-d2-ul2', FVU, 'pure subroutine gideon_fv_d2_ul2', 'd = (gl-gr)', 'd = (gr-gl)'),
    ('face-offset-fv-d1-ur3', FVU, 'pure subroutine gideon_fv_d1_ur3', 'gl = gl + FV_UFACE_S3(m)*q(m-1)',
     'gl = gl + FV_UFACE_S3(m)*q(m)'),
    # dispatch: a wrong kernel for an order, a number instead of NaN
    ('dispatch-wrong-kernel', DIS, 'pure subroutine gideon_fd_d1_c(', 'case(4) ; call gideon_fd_d1_c4',
     'case(4) ; call gideon_fd_d2_c4'),
    ('dispatch-no-nan', DIS, 'pure subroutine gideon_fv_face_ur(', 'case default ; f = GIDEON_NAN',
     'case default ; f = 0._R_P'),
    ('nan-constant', 'src/lib/gideon_kind_parameters.F90', 'GIDEON_NAN = transfer(9221120237041090560_I8P',
     '9221120237041090560_I8P', '4607182418800017408_I8P'),
    # 3-D operators: a curl component, a step, a Laplacian direction
    ('curl-component', O3D, 'pure subroutine gideon_fd_curl_c6', 'c(2) = dvx_dz - dvz_dx', 'c(2) = dvz_dx - dvx_dz'),
    ('divergence-step', O3D, 'pure subroutine gideon_fd_divergence_c4', 'call gideon_fd_d1_c4(vy_y, dxyz(2), dy)',
     'call gideon_fd_d1_c4(vy_y, dxyz(1), dy)'),
    ('laplacian-direction', O3D, 'pure subroutine gideon_fv_laplacian_c6', 'call gideon_fv_d2_c6(qz, dxyz(3), dz)',
     'call gideon_fv_d2_c6(qy, dxyz(3), dz)'),
    # phase 2: FD face derivatives, one-sided and biased FD face values (each mirror mutant reads inside its window)
    ('row-fd-face-d1-c4', T, 'FD_FACE_D1_S2(1:2)', '9._R_P/8._R_P', '1.1251_R_P'),
    ('row-fd-uface-3', T, 'FD_UFACE_S3(1:3)', '-5._R_P/4._R_P', '-1.2501_R_P'),
    ('row-fd-bface-5-tail', T, 'FD_BFACE_S2(-1:3)', '-5._R_P/128._R_P', '-5._R_P/127._R_P'),
    ('pair-sign-fd-face-d2-c4', FDC, 'pure subroutine gideon_fd_face_d2_c4', '(q(m) + q(1-m))', '(q(m) - q(1-m))'),
    ('scale-fd-face-d2-c4', FDC, 'pure subroutine gideon_fd_face_d2_c4', 'd = d/ds/ds', 'd = d/ds'),
    ('mirror-fd-face-ul3', FDU, 'pure subroutine gideon_fd_face_ul3', 'q(1-m)', 'q(m-3)'),
    ('mirror-fd-face-bl5', FDU, 'pure subroutine gideon_fd_face_bl5', 'q(1-m)', 'q(m-1)'),
    ('dispatch-wrong-kernel-bl', DIS, 'pure subroutine gideon_fd_face_bl(', 'case(5) ; call gideon_fd_face_bl5',
     'case(5) ; call gideon_fd_face_bl3'),
    # P2 (issue #2): the explicit families of phases 1 and 2.1
    ('row-face-diff-5', T, 'FACE_DIFF_S3(1:3)', '5._R_P', '5.001_R_P'),
    ('bdiss-sign-fd-br5', DSS, 'pure subroutine gideon_fd_face_diss_br5', '-FD_BDISS_S2*d', 'FD_BDISS_S2*d'),
    ('bdiss-constant-fv-3', T, 'FV_BDISS_S1 = ', '1._R_P/12._R_P', '1._R_P/11._R_P'),
    ('row-fv-bface-5', T, 'FV_BFACE_S2(-1:3)', '47._R_P/60._R_P', '0.7834_R_P'),
    ('mirror-fv-face-bl5', FVU, 'pure subroutine gideon_fv_face_bl5', 'q(1-m)', 'q(m-1)'),
    ('row-filter-c6', T, 'FILTER_S3(0:3)', '15._R_P/64._R_P', '0.2344_R_P'),
    ('pair-sign-fv-filter-c4', FIL, 'pure subroutine gideon_fv_filter_c4', '(q(m) + q(-m))', '(q(m) - q(-m))'),
    ('dispatch-wrong-kernel-diss', DIS, 'pure subroutine gideon_fd_face_diss_c(', 'case(6) ; call gideon_fd_face_diss_c6',
     'case(6) ; call gideon_fd_face_diss_c4'),
    # P2: the line solvers (the periodic ones through their bounded factor, the masked through their merges)
    ('thomas-pivot', LS, 'pure subroutine tri_factor_const', 'p = b - a*w(i-1)', 'p = b + a*w(i-1)'),
    ('tri-periodic-column', LS, 'pure subroutine gideon_tri_periodic_var', 'z(n-1) = z(n-1) + c(n-1)',
     'z(n-1) = z(n-1) + a(n-1)'),
    ('penta-resolve-pivot', LS, 'pure subroutine penta_resolve_const', 'bp = bp - ap*u1(i-1)', 'bp = bp + ap*u1(i-1)'),
    ('penta-border', LS, 'pure subroutine gideon_penta_periodic_const', 'x1 = (r1 - s12*x2)/s11',
     'x1 = (r1 + s12*x2)/s11'),
    ('masked-diagonal', LS, 'pure subroutine tri_factor_masked', 'p = merge(1._R_P, b, mask(i))', 'p = b'),
    ('masked-offdiagonal', LS, 'pure subroutine penta_factor_masked', 'ap = ap - merge(0._R_P, e, mask(i))',
     'ap = ap - e'),
    # P2: the compact operators (bands, right sides, signs, steps, the filters' alpha_f, the biased mirror)
    ('band-fd-d1-p8', T, 'FD_D1_P8_LEFT(0:2)', '4._R_P/9._R_P', '0.4445_R_P'),
    ('right-fd-d2-p6', T, 'FD_D2_P6_RIGHT(0:2)', '12._R_P/11._R_P', '1.0910_R_P'),
    ('pair-sign-fd-face-d1-p6', CMP, 'pure subroutine gideon_fd_face_d1_p6', '(q(m) - q(1-m))', '(q(m) + q(1-m))'),
    ('scale-fd-d2-p4', CMP, 'pure subroutine gideon_fd_d2_p4', 'r = r/ds/ds', 'r = r/ds'),
    ('filter-alpha-row-p6', T, 'FD_FILTER_P6_ALPHA(0:3)', '17._R_P/32._R_P', '0.5313_R_P'),
    ('filter-alpha-term-p4', CMP, 'pure subroutine gideon_fd_filter_p4', 'alpha*FD_FILTER_P4_ALPHA(m)',
     'alpha*FD_FILTER_P4_ALPHA(0)'),
    ('mirror-fv-face-pbl7', CMP, 'pure subroutine gideon_fv_face_pbl7', 'q(1-m)', 'q(m-1)'),
    ('band-fv-face-pbl5', T, 'FV_FACE_PBL5_LEFT(-1:1)', '1._R_P/2._R_P', '0.5001_R_P'),
    ('band-fv-face-p6', T, 'FV_FACE_P6_LEFT(0:1)', '1._R_P/3._R_P', '0.3334_R_P'),
    # #2 decision 4a: factor once, solve many (a pivot, the carry of row i-2, the periodic border)
    ('factor-pivot-tri', LS, 'pure subroutine tri_factor_rows', 'fac(i,1) = 1._R_P/(b - a*wp)',
     'fac(i,1) = 1._R_P/(b + a*wp)'),
    ('factor-carry-penta', LS, 'pure subroutine penta_factor_rows', 'u1pp = u1p ; u2pp = u2p',
     'u1pp = u1 ; u2pp = u2'),
    ('factored-border-penta', LS, 'pure subroutine gideon_penta_periodic_solve', 'x1 = (r1 - fac(n-1,2)*x2)*fac(n-1,1)',
     'x1 = (r1 + fac(n-1,2)*x2)*fac(n-1,1)'),
]


def run(cmd: list[str], cwd: Path, timeout: int = 3600) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)


def failing_tests(work: Path, mode: str) -> tuple[bool, list[str]]:
    """Build the tests; return (built, names of the failing test programs)."""
    # --force-compile: FoBiS does not rebuild the dependents of a changed cpp #include (the kernels' .INC bodies), so
    # an incremental build would keep a stale mutant compiled in, or miss a new one
    b = run(['fobis', 'build', '--mode', mode, '--force-compile'], work)
    if b.returncode != 0 or 'Error' in b.stdout and 'error' in b.stdout.lower() and not (work / 'exe').exists():
        return False, []
    t = run(['bash', 'scripts/run_tests.sh'], work)
    # the report lines of run_tests.sh: '  FAIL  <name> (<reason>)' (a test's own 'FAILED: ...' messages are not)
    failed = [line.split()[1].split('/')[-1] for line in t.stdout.splitlines()
              if line.split()[:1] in (['FAIL'], ['XPASS'])]
    if t.returncode != 0 and not failed:
        failed = ['(run_tests.sh failed)']
    return True, failed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--mode', default='tests-gnu')
    ap.add_argument('--workdir', type=Path)
    ap.add_argument('--only', nargs='*')
    args = ap.parse_args()
    repo = Path.cwd()
    work = args.workdir or Path(tempfile.mkdtemp(prefix='gideon-mutants-'))
    if work.exists():
        shutil.rmtree(work)
    files = run(['git', 'ls-files'], repo).stdout.split() + run(['git', 'ls-files', '--others', '--exclude-standard'],
                                                                repo).stdout.split()
    for f in files:
        if not (repo / f).exists():  # tracked but deleted in the working tree (a page moved, not yet committed)
            continue
        dst = work / f
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / f, dst)
    shutil.copytree(repo / 'src/third_party/PENF', work / 'src/third_party/PENF')
    built, failed = failing_tests(work, args.mode)
    if not built or failed:
        print(f'baseline is not green ({failed}): fix the suite first')
        return 1
    print(f'baseline green in {work}')
    survivors = []
    for name, rel, anchor, old, new in MUTANTS:
        if args.only and name not in args.only:
            continue
        path = work / rel
        original = path.read_text()
        i = original.find(anchor)
        j = original.find(old, i) if i >= 0 else -1
        if i < 0 or j < 0:
            print(f'{name:28s} NOT APPLIED (anchor or text not found)')
            survivors.append(name)
            continue
        path.write_text(original[:j] + new + original[j + len(old):])
        built, failed = failing_tests(work, args.mode)
        path.write_text(original)
        if not built:
            print(f'{name:28s} does not build (invalid mutant)')
            survivors.append(name)
        elif failed:
            print(f'{name:28s} killed by {", ".join(sorted(set(failed)))}')
        else:
            print(f'{name:28s} SURVIVED')
            survivors.append(name)
        sys.stdout.flush()
    print(f'{len(MUTANTS) - len(survivors)} of {len(MUTANTS)} mutants killed' if not args.only else 'done')
    if survivors:
        print('survivors: ' + ', '.join(survivors))
    return 0 if not survivors else 1


if __name__ == '__main__':
    sys.exit(main())
