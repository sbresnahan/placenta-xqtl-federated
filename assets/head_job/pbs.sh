#!/bin/bash
# Run the Nextflow head as a PBS Pro / Torque job.  Edit, then:  qsub pbs.sh
#PBS -N xqtl_head
#PBS -q long
#PBS -l walltime=120:00:00
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -j oe
set -euo pipefail
export NXF_OPTS="-Xms1g -Xmx4g"
cd "${PBS_O_WORKDIR:-$PWD}"
./xqtl/run_xqtl.sh -profile "${PROFILE:-pbspro,apptainer}" -c "${SITE_CONFIG:-xqtl/pipeline/conf/sites/mysite.config}" \
    -params-file "${COHORT_YML:-cohort.yml}" --mode "${MODE:-site}" -resume
