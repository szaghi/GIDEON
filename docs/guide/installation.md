# Installation

## FoBiS

```bash
git clone https://github.com/szaghi/GIDEON && cd GIDEON
fobis fetch                      # PENF, into src/third_party/
fobis build --mode static-gnu    # lib/gnu/static/libgideon.a (shared-gnu: libgideon.so)
```

`fobis build --lmodes` lists the modes: static and shared libraries and tests for gfortran (`-gnu`), nvfortran
(`-nvf`) and ifx (`-intel`), `tests-nvf-oac`, the OpenACC build (`--varset cc80|cc89|cc90` picks the compute
capability; cc89 by default), and the single-precision twins (`-r4`, `R_P = R4P`).

## fpm

```toml
[dependencies]
GIDEON = { git = "https://github.com/szaghi/GIDEON", branch = "main" }
```
