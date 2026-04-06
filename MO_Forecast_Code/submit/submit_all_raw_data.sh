#!/bin/bash
set -euo pipefail

SPEC_DIR="specs/raw_data"
PBS_SCRIPT="submit_one_raw_spec.pbs"

extract_ppn () {
  local yml="$1"
  local w

  w="$(grep -E '^[[:space:]]*workers:[[:space:]]*[0-9]+' -m 1 "$yml" \
      | sed -E 's/^[[:space:]]*workers:[[:space:]]*([0-9]+).*/\1/')"

  if [[ -z "${w}" ]]; then
    echo 26
  else
    echo $((w + 1))
  fi
}

for yml in "${SPEC_DIR}"/*.yml; do
  [[ -e "$yml" ]] || { echo "No YAML files found"; exit 1; }

  spec_id="$(basename "$yml" .yml)"
  ppn="$(extract_ppn "$yml")"

  echo "Submitting ${spec_id} with ppn=${ppn}"

  qsub \
    -N "onset_${spec_id}" \
    -v SPEC_ID="${spec_id}" \
    -l nodes=1:ppn="${ppn}" \
    -o cluster_outputs/ \
    -j oe \
    "${PBS_SCRIPT}"
done
