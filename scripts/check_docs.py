#!/usr/bin/env python3
"""Check that the documentation covers the code, and that it names nothing the code does not have.

Usage: check_docs.py [--verbose]

The pages are the hand-written ones under docs/ (docs/api, written by `formal`, and docs/node_modules excluded), with
the tables they include. Every check prints its count; a failure lists what is missing and the script exits 1.

Coverage (the code -> the pages):
  symbols    every public procedure, kernel, dispatcher, 3-D operator, row and kind of src/lib is named on a page
  modules    every module of src/lib is named on a page
  tests      every test program of src/tests is named on a page of the verification or regression sections
  examples   every program of docs/examples/src is included by a page, and its output too
  modes      every fobos mode is named on a page
  scripts    every script of scripts/ and docs/scripts/ is named on a page
  figures    every file of docs/public is referenced by a page (a -dark twin through its light figure)
  pages      every page is in the sidebar or linked from a page; every sidebar entry is a page

Accuracy (the pages -> the code):
  names      every gideon_* name on a page is a procedure, module or program, or one of their families written with a
             placeholder (dN, <p>, {ur,ul})
  adam       every ADAM routine named on a page exists in ADAM's verbatim copy (src/tests/adam_*.F90)
  links      every image and every internal link of a page resolves to a file
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / 'docs'
PUBLIC = DOCS / 'public'

FAILURES: list[str] = []


def report(check: str, total: int, missing: list[str], verbose: bool) -> None:
    """Print one check's result and record its failures."""
    status = 'ok' if not missing else f'FAIL ({len(missing)})'
    print(f'{check:9s} {total:5d} checked  {status}')
    for item in sorted(set(missing)):
        FAILURES.append(f'{check}: {item}')
        print(f'          missing: {item}')
    if verbose and not missing:
        print('          (all covered)')


def pages() -> list[Path]:
    """Every hand-written page and included table (not docs/ford.md, the API generator's project file)."""
    skip = {'api', 'node_modules', '.vitepress', 'public', 'examples'}
    return sorted(p for p in DOCS.rglob('*.md') if not skip & set(p.relative_to(DOCS).parts) and p.name != 'ford.md')


def text_of(paths: list[Path]) -> str:
    return '\n'.join(p.read_text(encoding='utf-8') for p in paths)


def fortran_sources(directory: Path) -> list[Path]:
    return sorted(directory.glob('*.F90'))


def library_symbols() -> dict[str, set[str]]:
    """The public entities of src/lib, by kind."""
    procs, modules, rows = set(), set(), set()
    for f in fortran_sources(ROOT / 'src' / 'lib'):
        src = f.read_text(encoding='utf-8')
        modules |= set(re.findall(r'^\s*module\s+(\w+)\s*$', src, re.M | re.I))
        procs |= set(re.findall(r'^\s*(?:pure\s+)?subroutine\s+(\w+)', src, re.M | re.I))
        procs |= set(re.findall(r'^\s*(?:pure\s+)?(?:elemental\s+)?function\s+(\w+)', src, re.M | re.I))
        if f.name == 'gideon_tables.F90':
            rows |= set(re.findall(r'parameter\s*::\s*(F[DV]_\w+)', src))
    kinds = {'R_P', 'GIDEON_NAN'}
    return {'procedures': procs, 'modules': modules, 'rows': rows, 'kinds': kinds}


def programs(directory: Path) -> set[str]:
    names = set()
    for f in fortran_sources(directory):
        names |= set(re.findall(r'^\s*program\s+(\w+)', f.read_text(encoding='utf-8'), re.M | re.I))
    return names


def fobos_modes() -> list[str]:
    text = (ROOT / 'fobos').read_text(encoding='utf-8')
    block = re.search(r'^\[modes\]\n(.*?)(?=^\[)', text, re.M | re.S)
    if block is None:
        return []
    body = re.sub(r'#.*', '', block.group(1)).replace('modes =', ' ')
    return body.split()


def sidebar_links() -> set[str]:
    config = (DOCS / '.vitepress' / 'config.mts').read_text(encoding='utf-8')
    return set(re.findall(r"link:\s*'(/[^']*)'", config))


def page_of_link(link: str) -> Path | None:
    """The page a site link resolves to (without its hash), or None."""
    path = link.split('#', 1)[0].split('?', 1)[0]
    if not path.startswith('/'):
        return None
    candidates = [DOCS / (path.strip('/') + '.md'), DOCS / path.strip('/') / 'index.md']
    if path == '/':
        candidates = [DOCS / 'index.md']
    for c in candidates:
        if c.is_file():
            return c
    return None


def check_symbols(text: str, verbose: bool) -> None:
    sym = library_symbols()
    words = set(re.findall(r'\b\w+\b', text))
    for kind in ('procedures', 'rows', 'kinds'):
        names = sym[kind]
        report(f'{kind[:9]}', len(names), [n for n in names if n not in words], verbose)
    report('modules', len(sym['modules']), [m for m in sym['modules'] if m not in words], verbose)


def check_tests(verbose: bool) -> None:
    section = [p for p in pages() if p.relative_to(DOCS).parts[0] in ('verification', 'regression')]
    words = set(re.findall(r'\b\w+\b', text_of(section)))
    names = programs(ROOT / 'src' / 'tests')
    report('tests', len(names), [n for n in names if n not in words], verbose)


def check_examples(text: str, verbose: bool) -> None:
    missing = []
    srcs = sorted(p for p in (DOCS / 'examples' / 'src').iterdir() if p.suffix in ('.F90', '.py'))
    for s in srcs:
        if f'examples/src/{s.name}' not in text:
            missing.append(f'docs/examples/src/{s.name} (no page includes it)')
    for out in sorted((DOCS / 'examples' / 'output').glob('*.txt')):
        if f'examples/output/{out.name}' not in text:
            missing.append(f'docs/examples/output/{out.name} (no page includes it)')
    report('examples', len(srcs), missing, verbose)


def check_modes(text: str, verbose: bool) -> None:
    modes = fobos_modes()
    report('modes', len(modes), [m for m in modes if not re.search(rf'(?<![\w-]){re.escape(m)}(?![\w-])', text)
                                 and not mode_family_named(m, text)], verbose)


def mode_family_named(mode: str, text: str) -> bool:
    """A mode named through a family, e.g. `static-{gnu,nvf,intel}-r4` or `shared-*-r4`."""
    for fam in re.findall(r'`([\w-]*\*[\w-]*)`', text):
        if re.fullmatch(re.escape(fam).replace(r'\*', r'[\w]+'), mode):
            return True
    for fam in re.findall(r'[\w-]*\{[\w,-]+\}[\w-]*', text):
        head, alts, tail = re.match(r'([\w-]*)\{([\w,-]+)\}([\w-]*)', fam).groups()
        if any(mode == f'{head}{a}{tail}' for a in alts.split(',')):
            return True
    return False


def check_scripts(text: str, verbose: bool) -> None:
    scripts = [p for p in (ROOT / 'scripts').iterdir() if p.is_file() and p.suffix in ('.py', '.sh')]
    scripts += [p for p in (ROOT / 'scripts' / 'ci').glob('*') if p.is_file()]
    scripts += [p for p in (DOCS / 'scripts').glob('*.py')]
    missing = [str(p.relative_to(ROOT)) for p in scripts if p.name not in text]
    report('scripts', len(scripts), missing, verbose)


def check_figures(text: str, verbose: bool) -> None:
    files = [p for p in PUBLIC.rglob('*') if p.is_file() and p.name != 'coverage.json']
    missing = []
    for f in files:
        rel = '/' + f.relative_to(PUBLIC).as_posix()
        light = re.sub(r'-dark\.png$', '.png', rel)
        if rel not in text and light not in text:
            missing.append(rel)
    report('figures', len(files), missing, verbose)


def check_pages(verbose: bool) -> None:
    all_pages = [p for p in pages() if not p.name.startswith('_') and '_tables' not in p.parts]
    sidebar = sidebar_links()
    reached = {page_of_link(link) for link in sidebar} - {None}
    linked = set()
    for p in all_pages:
        for link in re.findall(r'\]\((/[^)\s]*)\)', p.read_text(encoding='utf-8')):
            linked.add(page_of_link(link))
    missing = [str(p.relative_to(ROOT)) for p in all_pages if p not in reached | linked and p.name != 'index.md'
               or (p.name == 'index.md' and p not in reached | linked and p.parent != DOCS)]
    missing += [f'sidebar link {link} (no page)' for link in sidebar if page_of_link(link) is None
                and not link.startswith('/api')]
    report('pages', len(all_pages), missing, verbose)


def family_regex(name: str) -> re.Pattern[str]:
    """A name written with placeholders, as a pattern matching the start of a real name.

    `<p>`, `<data>`: any word; `{a,b}`: one of them; `{2..5}`: a number; `dN`, `derivativeN`: a derivative order. A
    family name matches the real names it begins (`gideon_fv_face_ur` the kernels `gideon_fv_face_ur1` ... and the
    dispatcher)."""
    out, i = [], 0
    for m in re.finditer(r'<\w+>|\{\d+\.\.\d+\}|\{[\w,]+\}|(?<=d)N(?![a-z])|(?<=derivative)N', name):
        out.append(re.escape(name[i:m.start()]))
        tok = m.group(0)
        if tok.startswith('<'):
            out.append(r'\w+')
        elif '..' in tok or tok == 'N':
            out.append(r'\d+')
        else:
            out.append('(' + '|'.join(re.escape(a) for a in tok[1:-1].split(',')) + ')')
        i = m.end()
    out.append(re.escape(name[i:]))
    # a family name begins its members; a plain name is itself or a dispatcher-like stem of kernels (`..._c` -> `..._c2`)
    tail = '' if i else r'(\d+)?$'
    if i == 0 and not re.search(r'_(c|ur|ul)$', name):
        tail = '$'
    return re.compile('^' + ''.join(out) + tail)


def check_names(text: str, verbose: bool) -> None:
    known = (library_symbols()['procedures'] | library_symbols()['modules'] | programs(ROOT / 'src' / 'tests')
             | programs(ROOT / 'src' / 'bench') | programs(DOCS / 'examples' / 'src'))
    known |= {f.stem for f in fortran_sources(ROOT / 'src' / 'tests')}
    known |= {'gideon'} | {f.stem for f in (DOCS / 'scripts').glob('*.py')}
    found = set(re.findall(r'gideon_(?:[\w<>]|\{[\w,.]+\})+[…*]?', text))
    missing = []
    for n in sorted(found):
        if n in known:
            continue
        if n[-1] in '…*':  # a prefix: gideon_test_…, gideon_test_*
            stem = family_regex(n[:-1]).pattern.rstrip('$').removesuffix(r'(\d+)?')
            if any(re.match(stem, k) for k in known):
                continue
            missing.append(n)
            continue
        if any(family_regex(n).match(k) for k in known):
            continue
        if re.match(r'gideon_golden_(r4_|r16_)?v\d+', n) or n in {'gideon_r4', 'gideon_r16'}:
            continue
        missing.append(n)
    report('names', len(found), missing, verbose)


def check_adam(text: str, verbose: bool) -> None:
    adam = ''.join(f.read_text(encoding='utf-8') for f in (ROOT / 'src' / 'tests').glob('adam_*.F90'))
    routines = set(re.findall(r'^\s*(?:pure\s+)?subroutine\s+(\w+)', adam, re.M | re.I))
    named = set(re.findall(r'\bcompute_(?:\w|\{[\w,.]+\})+', text))
    missing = [n for n in named if n not in routines and not any(family_regex(n).match(r) for r in routines)]
    report('adam', len(named), missing, verbose)


def check_links(verbose: bool) -> None:
    missing, total = [], 0
    for p in pages():
        src = p.read_text(encoding='utf-8')
        for img in re.findall(r'!\[[^\]]*\]\((/[^)\s]+)\)', src):
            total += 1
            if not (PUBLIC / img.lstrip('/')).is_file():
                missing.append(f'{p.relative_to(ROOT)}: image {img}')
        for media in re.findall(r'(?:src|poster)="(/[^"]+)"', src):
            total += 1
            if not (PUBLIC / media.lstrip('/')).is_file():
                missing.append(f'{p.relative_to(ROOT)}: media {media}')
        for link in re.findall(r'(?<!!)\[[^\]]*\]\((/[^)\s]*)\)', src):
            total += 1
            if link.startswith('/api') or (PUBLIC / link.lstrip('/').split('#')[0]).is_file():
                continue
            if page_of_link(link) is None:
                missing.append(f'{p.relative_to(ROOT)}: link {link}')
        for inc in re.findall(r'<!--@include:\s*(\S+?)\s*-->', src):
            total += 1
            if not (p.parent / inc).resolve().is_file():
                missing.append(f'{p.relative_to(ROOT)}: include {inc}')
        for snip in re.findall(r'<<<\s*@/(\S+?)(?:\{|\s|$)', src):
            total += 1
            if not (DOCS / snip).is_file():
                missing.append(f'{p.relative_to(ROOT)}: snippet {snip}')
    report('links', total, missing, verbose)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split('\n', 1)[0])
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()
    text = text_of(pages())
    check_symbols(text, args.verbose)
    check_tests(args.verbose)
    check_examples(text, args.verbose)
    check_modes(text, args.verbose)
    check_scripts(text, args.verbose)
    check_figures(text, args.verbose)
    check_pages(args.verbose)
    check_names(text, args.verbose)
    check_adam(text, args.verbose)
    check_links(args.verbose)
    if FAILURES:
        print(f'check_docs: {len(FAILURES)} failure(s)')
        return 1
    print('check_docs: ok')
    return 0


if __name__ == '__main__':
    sys.exit(main())
