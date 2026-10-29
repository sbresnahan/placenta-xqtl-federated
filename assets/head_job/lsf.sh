#!/bin/bash
# Run the Nextflow head as an LSF job (docs/HPC.md, "Running the head inside a batch job").
# Edit queue/walltime/profile, then:  bsub < lsf.sh
#BSUB -J xqtl_head
#BSUB -q long
#BSUB -W 120:00
#BSUB -n 1
#BSUB -M 8
#BSUB -R "rusage[mem=8]"
#BSUB -o xqtl_head.%J.out
#BSUB -e xqtl_head.%J.err
set -euo pipefail
export NXF_OPTS="-Xms1g -Xmx4g"
cd "${LS_SUBCWD:-$PWD}"
MODE=${MODE:-site}
./xqtl/run_xqtl.sh -profile "${PROFILE:-seadragon}" -params-file "${COHORT_YML:-cohort.yml}" --mode "$MODE" -resume
