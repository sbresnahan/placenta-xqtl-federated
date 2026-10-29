#!/bin/bash
# Run the Nextflow head as a SLURM job.  Edit, then:  sbatch slurm.sh
#SBATCH --job-name=xqtl_head
#SBATCH --partition=long
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=xqtl_head.%j.out
set -euo pipefail
export NXF_OPTS="-Xms1g -Xmx4g"
cd "${SLURM_SUBMIT_DIR:-$PWD}"
./xqtl/run_xqtl.sh -profile "${PROFILE:-slurm,apptainer}" -c "${SITE_CONFIG:-xqtl/pipeline/conf/sites/mysite.config}" \
    -params-file "${COHORT_YML:-cohort.yml}" --mode "${MODE:-site}" -resume
