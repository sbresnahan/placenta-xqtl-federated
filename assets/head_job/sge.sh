#!/bin/bash
# Run the Nextflow head as an SGE/UGE job.  Edit, then:  qsub sge.sh
#$ -N xqtl_head
#$ -cwd
#$ -q long.q
#$ -l h_rt=120:00:00
#$ -l h_vmem=8G
#$ -j y
set -euo pipefail
export NXF_OPTS="-Xms1g -Xmx4g"
./xqtl/run_xqtl.sh -profile "${PROFILE:-sge,apptainer}" -c "${SITE_CONFIG:-xqtl/pipeline/conf/sites/mysite.config}" \
    -params-file "${COHORT_YML:-cohort.yml}" --mode "${MODE:-site}" -resume
