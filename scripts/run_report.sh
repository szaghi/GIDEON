#!/usr/bin/env bash
# The data, tables, figures and animations of the verification and regression pages, from GIDEON's own programs.
#
# Usage: scripts/run_report.sh [--no-build] [--gpu] [--report-dir DIR] [--skip-animations]
#   --no-build         use the test programs already built (exe/, exe-r16/, exe-r4/)
#   --gpu              also build the OpenACC tests (tests-nvf-oac into exe/oac/) and run gideon_test_device on the GPU
#                      of CUDA_VISIBLE_DEVICES (default 1); needs nvfortran on PATH
#   --report-dir DIR   where the data go (default: exe/report); a mutation run saved there as mutation-run-<date>.txt
#                      (python3 scripts/mutation_test.py > DIR/mutation-run-<date>.txt) is tabled too
#   --skip-animations  do not render the casts (ffmpeg)
#
# Steps:
#   1. build the test programs: fobis build --mode tests-gnu, tests-gnu-r16 (the orders observed in quadruple
#      precision), tests-gnu-r4 (the single-precision golden values);
#   2. the generator's rows, windows, orders and findings: scripts/gen_gideon.py --check --json rows.json;
#   3. every test program of each build, with GIDEON_REPORT_DIR set: the programs write their data (CSV, frames) and
#      their output is kept (<build>/<program>.txt); a failing program stops the report;
#   4. docs/scripts/make_tables.py (docs/verification/_tables, docs/regression/_tables), make_figures.py
#      (docs/public/verification, docs/public/guide, light and dark), make_animations.py (docs/public/verification/casts).
#
# Python: the interpreter of GIDEON_DOCS_PYTHON, else exe/.docs-venv/bin/python when it exists, else python3; it needs
# numpy and matplotlib (and ffmpeg on PATH for the animations).
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

BUILD=1 ; GPU=0 ; ANIM=1
REPORT="$ROOT/exe/report"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-build ) BUILD=0; shift ;;
    --gpu ) GPU=1; shift ;;
    --report-dir ) REPORT="$2"; shift 2 ;;
    --skip-animations ) ANIM=0; shift ;;
    * ) printf "Unknown argument: %s\n" "$1" >&2; exit 2 ;;
  esac
done
if [[ -n "${GIDEON_DOCS_PYTHON:-}" ]]; then
  PY="$GIDEON_DOCS_PYTHON"
elif [[ -x "$ROOT/exe/.docs-venv/bin/python" ]]; then
  PY="$ROOT/exe/.docs-venv/bin/python"
else
  PY=python3
fi

if [[ $BUILD -eq 1 ]]; then
  fobis build --mode tests-gnu
  fobis build --mode tests-gnu-r16
  fobis build --mode tests-gnu-r4
  [[ $GPU -eq 0 ]] || fobis build --mode tests-nvf-oac --build-dir exe/oac/
fi

mkdir -p "$REPORT"
rm -rf "$REPORT/r32" "$REPORT/r64" "$REPORT/r128"
if [[ $GPU -eq 1 ]]; then
  rm -rf "$REPORT/oac"
elif [[ -d "$REPORT/oac" ]]; then
  printf "warning: %s/oac is the OpenACC data of an earlier run (%s): rerun with --gpu after a change of the kernels\n" \
    "$REPORT" "$(date -r "$REPORT/oac/gideon_test_device.txt" +%F 2>/dev/null || echo '?')" >&2
fi
python3 scripts/gen_gideon.py --check --json "$REPORT/rows.json" > "$REPORT/gen_gideon.txt"
cat "$REPORT/gen_gideon.txt"

# run every program of a build directory into $REPORT/<tag>
run_build() {
  local dir="$1" tag="$2" exe name
  mkdir -p "$REPORT/$tag"
  for exe in "$dir"/*; do
    [[ -f "$exe" && -x "$exe" ]] || continue
    name=$(basename "$exe")
    if ! GIDEON_REPORT_DIR="$REPORT/$tag" "$exe" > "$REPORT/$tag/$name.txt" 2>&1; then
      printf "FAIL %s/%s:\n" "$tag" "$name" >&2
      cat "$REPORT/$tag/$name.txt" >&2
      exit 1
    fi
    printf "  ran %-8s %s\n" "$tag" "$name"
  done
}
run_build exe r64
run_build exe-r16 r128
run_build exe-r4 r32
if [[ $GPU -eq 1 ]]; then
  export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"
  run_build exe/oac oac
  nvidia-smi --query-gpu=name,driver_version --format=csv,noheader -i "$CUDA_VISIBLE_DEVICES" \
    > "$REPORT/oac/gpu.txt" || true
  nvfortran --version | sed -n '2p' > "$REPORT/oac/compiler.txt" || true
  cat "$REPORT/oac/gideon_test_device.txt"
fi
gfortran --version | head -1 > "$REPORT/compiler.txt"

"$PY" docs/scripts/make_tables.py --report "$REPORT"
"$PY" docs/scripts/make_figures.py --report "$REPORT" --theme both
if [[ $ANIM -eq 1 ]]; then
  "$PY" docs/scripts/make_animations.py --report "$REPORT"
fi
