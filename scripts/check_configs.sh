#!/usr/bin/env bash
# Resolve the configuration for every supported profile combination and check key values.
# Fails if any profile does not parse or resolves to an unexpected executor/engine.
#   scripts/check_configs.sh            (run from anywhere; needs nextflow on PATH)
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
export NXF_OFFLINE=true
fail=0
check() {  # profile, expected executor, expected engine (docker|singularity|apptainer|none)
    local prof=$1 exe=$2 eng=$3 out
    if ! out=$(nextflow config -flat -profile "$prof" . 2>&1); then
        echo "FAIL  $prof : does not resolve"; echo "$out" | tail -5 | sed 's/^/      /'; fail=1; return
    fi
    local got_exe got_eng=none
    got_exe=$(echo "$out" | awk -F' = ' '$1=="process.executor"{gsub(/\x27/,"",$2); print $2}')
    for e in docker singularity apptainer; do
        echo "$out" | grep -q "^$e.enabled = true" && got_eng=$e
    done
    if [ "${got_exe:-local}" = "$exe" ] && [ "$got_eng" = "$eng" ]; then
        echo "OK    $prof  (executor=${got_exe:-local}, engine=$got_eng)"
    else
        echo "FAIL  $prof : executor=${got_exe:-local} (want $exe), engine=$got_eng (want $eng)"; fail=1
    fi
}
check seadragon                lsf    singularity
check seadragon,smoke          lsf    singularity
check lsf,apptainer            lsf    apptainer
check slurm,apptainer          slurm  apptainer
check slurm,singularity        slurm  singularity
check pbspro,apptainer         pbspro apptainer
check pbs,apptainer            pbs    apptainer
check sge,apptainer            sge    apptainer
check local,apptainer          local  apptainer
check local,dev,smoke          local  none
check test,docker              local  docker
check test,dev                 local  none
check test_phenotyping,docker  local  docker
check test_phenotyping,dev     local  none
for t in conf/sites/TEMPLATE_*.config; do
    sched=$(basename "$t" .config | sed 's/TEMPLATE_//')
    if out=$(nextflow -c "$t" config -flat -profile "$sched,apptainer" . 2>&1); then
        echo "OK    -c $t -profile $sched,apptainer"
        # values set in the site file must reach the container binds and the resource caps
        if echo "$out" | grep -q "^apptainer.runOptions = '-B /scratch -B /project'" \
           && echo "$out" | grep -q "^process.resourceLimits.memory = '256 GB'"; then
            echo "OK    -c $t : bind paths and resource caps from the site file"
        else
            echo "FAIL  -c $t : bind paths or resource caps from the site file not applied"; fail=1
        fi
    else
        echo "FAIL  -c $t"; echo "$out" | tail -5 | sed 's/^/      /'; fail=1
    fi
done
# seadragon: bind paths and resource caps from the site profile
sd=$(nextflow config -flat -profile seadragon .)
if echo "$sd" | grep -q "^singularity.runOptions = '-B /rsrch5 -B /rsrch9'" \
   && echo "$sd" | grep -q "^process.resourceLimits.memory = '120 GB'"; then
    echo "OK    seadragon bind paths and resource caps"
else
    echo "FAIL  seadragon bind paths or resource caps"; fail=1
fi
# Queue selection by walltime on seadragon (short <= 3h, medium < 24h, long <= 120h)
q=$(nextflow config -flat -profile seadragon . | grep -c "process.queue")
[ "$q" -ge 1 ] && echo "OK    seadragon queue rule present" || { echo "FAIL  seadragon queue rule"; fail=1; }
exit $fail
