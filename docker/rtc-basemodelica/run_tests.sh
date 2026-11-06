#!/bin/bash
# Run the RTC-Tools tests on all routes and print a summary.
#   run_tests.sh [unit] [examples] [examples-bm] [examples-fmu] [omc]    default: all of them
# Exit code: the number of failed parts.
set -u

RTC=${RTC_TOOLS:-/src/rtc-tools}
HERE=$(cd "$(dirname "$0")" && pwd)
EXAMPLE_TIMEOUT=${EXAMPLE_TIMEOUT:-1800}
cd "$RTC" || exit 1
parts=("$@")
[ ${#parts[@]} -eq 0 ] && parts=(unit examples examples-bm examples-fmu omc)
summary=()
failed=0

report() {
  summary+=("$1: $2")
  case $2 in
    *FAILED*) failed=$((failed + 1)) ;;
  esac
}

# run_scripts <part> <script>...: run each script from its own folder, as a user would
run_scripts() {
  local part=$1 passed=0 failures=() script started rc log
  shift
  log=$(mktemp)
  for script in "$@"; do
    started=$SECONDS
    (cd "$(dirname "$script")" && timeout "$EXAMPLE_TIMEOUT" python "$(basename "$script")") > "$log" 2>&1
    rc=$?
    echo "$part $script rc=$rc $((SECONDS - started))s"
    if [ $rc -eq 0 ]; then
      passed=$((passed + 1))
    else
      failures+=("$script")
      tail -3 "$log" | cut -c1-200
    fi
  done
  rm -f "$log"
  if [ ${#failures[@]} -eq 0 ]; then
    report "$part" "$passed passed"
  else
    report "$part" "$passed passed, FAILED: ${failures[*]}"
  fi
}

# The ModelicaMixin examples, which tests/examples/test_examples.py runs too
plain_examples() {
  ls examples/*/src/*example*.py | grep -v '_fmu\.py$'
}

for part in "${parts[@]}"; do
  echo "=== $part"
  case $part in
    unit)
      if python -m pytest --ignore=tests/examples --ignore=tests/modelica tests -q -p no:cacheprovider; then
        report unit passed
      else
        report unit FAILED
      fi
      ;;
    examples)
      # shellcheck disable=SC2046
      run_scripts examples $(plain_examples)
      ;;
    examples-bm)
      # shellcheck disable=SC2046
      PYTHONPATH="$HERE:${PYTHONPATH:-}" run_scripts examples-bm $(plain_examples)
      ;;
    examples-fmu)
      # shellcheck disable=SC2046
      run_scripts examples-fmu $(ls examples/*/src/*_fmu.py)
      ;;
    omc)
      # The ratchet compares against the omc version it declares: use the one in this image.
      version=$(omc --version | sed 's/^OpenModelica //')
      if OMC_VERSION=$version python -m pytest tests/modelica -q -p no:cacheprovider; then
        report omc passed
      else
        report omc FAILED
      fi
      ;;
    *)
      echo "unknown part: $part" >&2
      exit 2
      ;;
  esac
done

echo
echo "=== summary"
printf '%s\n' "${summary[@]}"
exit $failed
