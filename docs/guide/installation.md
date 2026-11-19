# Installation

GIDEON is a Fortran 2018 library with one dependency, [PENF](https://github.com/szaghi/PENF) (the kinds), compiled
from its sources by the same build. Its tests have two more: [SELAH](https://github.com/szaghi/SELAH), the source of
the time integrators a test needs (GIDEON implements space operators only), and [FUNDAL](https://github.com/szaghi/FUNDAL),
which SELAH uses. The library modes do not compile them. Three ways to build it: FoBiS (the project's build tool, every mode below), fpm,
or the installer script of a release.

## Toolchains

| Compiler | Version tested | Modes |
|---|---|---|
| gfortran | 16.2.0 | `*-gnu*`: libraries, tests (optimised, debug, coverage), R4P, R16P, benchmark, examples |
| nvfortran | 26.9 (NVIDIA HPC SDK) | `*-nvf*`: libraries, tests, the OpenACC tests and benchmark |
| ifx | 2026.2 (locally, and the CI matrix's `linux-oneapi-2026.2` image) | `*-intel*`: libraries, tests |

CI (`.github/workflows/`) builds and tests every push: `ci.yml` runs the coverage rule with gfortran (R8P, then R16P),
`matrix.yml` the tests of each compiler family (`tests-gnu`, `tests-nvf`, `tests-intel`) in the
[fortran-toolchains](https://github.com/szaghi/fortran-toolchains) images, `fpm.yml` the fpm build and tests;
`install.yml` installs each published release as its users do. The OpenACC modes need an NVIDIA GPU and are run by
hand.

## FoBiS

```bash
git clone https://github.com/szaghi/GIDEON && cd GIDEON
fobis fetch                        # PENF, and for the tests FUNDAL and SELAH, into src/third_party/
fobis build                        # the first mode, static-gnu: lib/gnu/static/libgideon.a
fobis build --mode tests-gnu && bash scripts/run_tests.sh
```

A program then compiles against the modules and links the archive, e.g.
`gfortran -I lib/gnu/static/mod prog.F90 lib/gnu/static/libgideon.a` (PENF's objects are inside the archive).
`fobis build --lmodes` lists every mode with its help line.

### Library modes

| Mode | Output | Flags (beyond the compiler's FoBiS profile) |
|---|---|---|
| `static-gnu`, `shared-gnu` | `lib/gnu/{static,shared}/libgideon.{a,so}` | `-std=f2018`, release profile `-O3 -funroll-loops`; shared: `-fPIC` |
| `static-nvf`, `shared-nvf` | `lib/nvf/{static,shared}/` | release profile `-O4 -fast`; shared: `-fPIC` |
| `static-intel`, `shared-intel` | `lib/intel/{static,shared}/` | `-O3 -standard-semantics -std18 -fp-model=precise`; shared: `-fPIC` |
| `static-*-r4`, `shared-*-r4` | `lib/<compiler>/{static,shared}-r4/libgideon_r4.{a,so}` | the same, and `-D_R_P_IS_R4P` |

The modules land in `mod/` under the output directory. The Intel modes do not use FoBiS's release profile: its `-ipo`
writes IR objects that FoBiS's link step and a plain `ar` archive cannot read, so `-O3` is given explicitly.

### Test modes

Every program of `src/tests/` into `exe/` (`exe-r4/`, `exe-r16/` for the other kinds), run by
`scripts/run_tests.sh [--exe-dir DIR]`:

| Mode | Purpose | Flags |
|---|---|---|
| `tests-gnu` | the default suite (CI matrix: family gnu) | `-std=f2018 -O3 -funroll-loops` |
| `tests-gnu-debug` | all run-time checks | debug profile: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra` |
| `tests-gnu-coverage` | coverage instrumentation (rule `makecoverage`) | coverage profile: `-O0 -g --coverage` |
| `tests-nvf` | the suite, nvfortran, CPU (CI matrix: family nvf) | `-O4 -fast -Kieee -tp=haswell` |
| `tests-nvf-debug` | all run-time checks | debug profile: `-O0 -g -C -traceback -Minfo=all`, `-Kieee -tp=haswell` |
| `tests-nvf-oac` | the suite with the device test on the GPU; `gideon_vv_advection`, a host program, links SELAH built for the host (`static-selah-nvf`) | `-acc -gpu=cc89 -Minfo=accel -Kieee -tp=haswell` |
| `tests-intel` | the suite, ifx (CI matrix: family intel) | `-O3 -standard-semantics -std18 -fp-model=precise` |
| `tests-intel-debug` | all run-time checks | debug profile: `-O0 -g -check all -traceback -warn all`, `-fp-model=precise` |
| `tests-gnu-r4`, `tests-gnu-debug-r4`, `tests-nvf-r4`, `tests-nvf-oac-r4`, `tests-intel-r4` | the same in single precision, into `exe-r4/` (without the ADAM comparison) | `-D_R_P_IS_R4P` |
| `static-selah-nvf`, `static-selah-nvf-r4` | [SELAH](https://szaghi.github.io/SELAH/)'s library from `src/third_party/SELAH`, nvfortran, host (no `-acc`), into `lib/nvf/`: built by the two OpenACC test modes (FoBiS `dependon`), whose host program `gideon_vv_advection` steps with it | `-Kieee -tp=haswell` |
| `tests-gnu-r16` | quadruple precision, host only, into `exe-r16/`: the orders of every kernel (without `gideon_vv_advection`: SELAH has no R16P build) | `-DPENF_R16P -D_R_P_IS_R16P` |

Every test mode adds FUNDAL's include directory (SELAH's kernels include `fundal.H`) and `-DFUNDAL_NO_MPI`. The two
OpenACC modes do not compile SELAH themselves: under `-acc` SELAH's directives need its device build, whose registers
are device memory, and `gideon_vv_advection`'s field is host memory; they link SELAH's host library instead.

Why `-Kieee -tp=haswell` and `-fp-model=precise`: the tests compare values bitwise, and both compilers otherwise
reorder the arithmetic, nvfortran differently on each CPU that compiles ([Single and quadruple precision](/cookbook/precision#the-flags-of-the-bitwise-checks)).

The OpenACC modes compile for compute capability 8.9 (an RTX 40xx) by default; `--varset cc80` (A100) or
`--varset cc90` (H100) picks another: `fobis build --mode tests-nvf-oac --varset cc90`.

### Other modes and rules

| Mode or rule | Purpose |
|---|---|
| `bench-nvf-oac`, `bench-nvf`, `bench-gnu` | `src/bench/gideon_bench_adam.F90` into `exe-bench/`: GIDEON's kernels against ADAM's on the GPU (`-acc -gpu=$GPU_CC`, `-fast` alone) and on one CPU core ([Performance](/guide/performance)) |
| `examples-gnu` | every program of `docs/examples/src/` into `exe-examples/`, checked against `docs/examples/output/` by `scripts/run_examples.sh` |
| `examples-gnu-r4` | the same programs in R4P into `exe-examples-r4/`, for the outputs `docs/examples/output/<name>-r4.txt` ([Precision](/cookbook/precision)); `scripts/run_examples.sh` builds it when such an output exists |
| `fobis rule --ex makecoverage` | clean, build `tests-gnu-coverage`, run the tests |
| `fobis rule --ex makecoverage-analysis` | the above, the coverage badge into `docs/public/coverage.json`, then `tests-gnu-r16` and its tests (what `ci.yml` runs) |
| `fobis rule --ex makedoc` | this site: `npm install` and `npm run docs:build` in `docs/` |
| `fobis rule --ex clean` | remove the build products and the generated documentation artefacts |

## fpm

```toml
[dependencies]
GIDEON = { git = "https://github.com/szaghi/GIDEON", branch = "main" }
```

The library is `src/lib`, in R8P; the tests are the `[[test]]` programs of `fpm.toml` (`fpm test`). The R4P and
R16P builds are FoBiS modes.

## The installer

Every release attaches `scripts/install.sh`: it downloads the release tarball (dependencies included) and builds and
installs it with the build system of your choice.

```bash
curl -fsSL https://github.com/szaghi/GIDEON/releases/latest/download/install.sh | bash -s -- --gnu --prefix ~/.local
```

| Option | Meaning |
|---|---|
| `--build fobis\|fpm\|cmake\|make` | the build system, among those the release ships |
| `--fc FC`, `--gnu`, `--intel`, `--nvidia` | the compiler (default `$FC`, else gfortran); the FoBiS modes are `static-gnu`, `static-intel`, `static-nvf` |
| `--prefix DIR` | where to install (default `~/.local`): `lib/gideon/<compiler-version>/`, `include/gideon/<compiler-version>/`, a pkg-config file |
| `--tag TAG`, `--source DIR` | a given release, or a local tree instead of a download |
| `--dry-run`, `--uninstall`, `--keep`, `--verbose` | check only; remove an install; keep the build directory; show the build output |

`bash scripts/install.sh --help` prints the full list.
