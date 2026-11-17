#!/usr/bin/env bash
# Build the documentation's example programs, run them and check their outputs.
#
# Usage: run_examples.sh [--write] [NAME ...]
#   Every program docs/examples/src/<name>.F90 is built by `fobis build --mode examples-gnu` into exe-examples/<name>,
#   run in a scratch directory, and its standard output compared with docs/examples/output/<name>.txt, the file the
#   tutorials and the cookbook include. A Python program docs/examples/src/<name>.py is run by `python3 -I` the same way. Exit 1 on a build failure, a non-zero exit status, a missing or different
#   output, or an output file without a program.
#   -w, --write  write docs/examples/output/<name>.txt from the runs instead of comparing (and remove the outputs of
#                programs that no longer exist).
#   NAME ...     only these programs (docs/examples/src/NAME.F90 or NAME.py); the check for stale outputs is skipped.
#
# A program with an output docs/examples/output/<name>-r4.txt is also built in single precision
# (`fobis build --mode examples-gnu-r4`, exe-examples-r4/) and that run compared with it; --write rewrites the R4P
# outputs that exist (create the empty file to opt a program in).
#
# The programs must print only what does not depend on the machine: explicit formats, no timings, no device names,
# no random numbers without a fixed seed.
set -uo pipefail

WRITE=0
NAMES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --write | -w ) WRITE=1; shift ;;
    -* ) printf "Unknown argument: %s\n" "$1" >&2; exit 2 ;;
    * ) NAMES+=("$1"); shift ;;
  esac
done

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
src_dir=$root/docs/examples/src
out_dir=$root/docs/examples/output
exe_dir=$root/exe-examples

if [[ -t 1 ]]; then
  RED=$'\033[0;31m'; GREEN=$'\033[0;32m'; BOLD=$'\033[1m'; RESET=$'\033[0m'
else
  RED=''; GREEN=''; BOLD=''; RESET=''
fi

scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT

cd "$root" || exit 1
[[ -d src/third_party/PENF/src ]] || fobis fetch > "$scratch/fetch.log" 2>&1 || {
  cat "$scratch/fetch.log"; echo "run_examples: fobis fetch failed" >&2; exit 1; }
fobis build --mode examples-gnu > "$scratch/build.log" 2>&1 || {
  cat "$scratch/build.log"; echo "run_examples: fobis build --mode examples-gnu failed" >&2; exit 1; }
shopt -s nullglob
r4_outputs=("$out_dir"/*-r4.txt)
if [[ ${#r4_outputs[@]} -gt 0 ]]; then
  fobis build --mode examples-gnu-r4 > "$scratch/build-r4.log" 2>&1 || {
    cat "$scratch/build-r4.log"; echo "run_examples: fobis build --mode examples-gnu-r4 failed" >&2; exit 1; }
fi

mkdir -p "$out_dir" "$scratch/run"
pass=0; fail=0
if [[ ${#NAMES[@]} -gt 0 ]]; then
  sources=()
  for name in "${NAMES[@]}"; do
    name=${name%.F90}; name=${name%.py}
    if [[ -f "$src_dir/$name.py" ]]; then sources+=("$src_dir/$name.py"); else sources+=("$src_dir/$name.F90"); fi
  done
else
  sources=("$src_dir"/*.F90 "$src_dir"/*.py)
fi
runs=()
for src in "${sources[@]}"; do
  name=$(basename "$src"); name=${name%.F90}; name=${name%.py}
  [[ $src == *.py ]] && { runs+=("$name:python:$name"); continue; }
  runs+=("$name:$exe_dir:$name")
  [[ -f "$out_dir/$name-r4.txt" ]] && runs+=("$name:$root/exe-examples-r4:$name-r4")
done
for run in "${runs[@]}"; do
  IFS=: read -r name exe_dir label <<< "$run"
  src=$src_dir/$name.F90
  [[ $exe_dir == python ]] && src=$src_dir/$name.py
  if [[ ! -f "$src" ]]; then
    printf "  ${RED}FAIL${RESET}  %s ${BOLD}(no program docs/examples/src/%s.F90)${RESET}\n" "$name" "$name"
    fail=$((fail + 1)); continue
  fi
  expected=$out_dir/$label.txt
  if [[ $exe_dir == python ]]; then
    cmd=(python3 -I "$src")
  else
    cmd=("$exe_dir/$name")
  fi
  if [[ $exe_dir != python && ! -x "$exe_dir/$name" ]]; then
    printf "  ${RED}FAIL${RESET}  %s ${BOLD}(not built: no %s)${RESET}\n" "$label" "${exe_dir#"$root"/}/$name"
    fail=$((fail + 1)); continue
  fi
  (cd "$scratch/run" && env LC_ALL=C OMP_NUM_THREADS=1 "${cmd[@]}" < /dev/null) > "$scratch/out" 2> "$scratch/err"
  rc=$?
  if [[ $rc -ne 0 ]]; then
    printf "  ${RED}FAIL${RESET}  %s ${BOLD}(exit status %d)${RESET}\n" "$label" "$rc"
    sed 's/^/       /' "$scratch/out" "$scratch/err"
    fail=$((fail + 1)); continue
  fi
  if [[ $WRITE -eq 1 ]]; then
    cp "$scratch/out" "$expected"
    printf "  ${GREEN}WROTE${RESET} %s\n" "${expected#"$root"/}"
    pass=$((pass + 1))
  elif [[ ! -f "$expected" ]]; then
    printf "  ${RED}FAIL${RESET}  %s ${BOLD}(no %s: run with --write)${RESET}\n" "$label" "${expected#"$root"/}"
    fail=$((fail + 1))
  elif cmp -s "$scratch/out" "$expected"; then
    printf "  ${GREEN}PASS${RESET}  %s\n" "$label"
    pass=$((pass + 1))
  else
    printf "  ${RED}FAIL${RESET}  %s ${BOLD}(output differs from %s)${RESET}\n" "$label" "${expected#"$root"/}"
    diff -u "$expected" "$scratch/out" | sed 's/^/       /'
    fail=$((fail + 1))
  fi
done

# an output without a program is stale: the page that includes it shows what no program prints
[[ ${#NAMES[@]} -gt 0 ]] && outputs=() || outputs=("$out_dir"/*.txt)
for out in "${outputs[@]}"; do
  name=$(basename "$out" .txt)
  name=${name%-r4}
  [[ -f "$src_dir/$name.F90" || -f "$src_dir/$name.py" ]] && continue
  if [[ $WRITE -eq 1 ]]; then
    rm -f -- "$out"
    printf "  ${GREEN}REMOVED${RESET} %s\n" "${out#"$root"/}"
  else
    printf "  ${RED}FAIL${RESET}  %s ${BOLD}(no program docs/examples/src/%s.F90)${RESET}\n" "${out#"$root"/}" "$name"
    fail=$((fail + 1))
  fi
done

total=$((pass + fail))
printf "\n${BOLD}%d/%d examples passed${RESET}\n" "$pass" "$total"
[[ $fail -eq 0 ]] || exit 1
