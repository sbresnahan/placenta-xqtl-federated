/*
 * --mode smoke : validate a new site before any real data are touched.
 *
 * Runs the Phase 0 checks (seadragon GO, 2026-10-29) against the single pipeline
 * image: container start on compute nodes, embedded tools and references, the
 * restricted file, CPU/memory enforcement, retry on memory kill, fan-out, job
 * arrays, bind mounts, tensorQTL (CPU, optional GPU) and the q-value bridge.
 * Every check writes "test_id <TAB> status <TAB> detail"; SMOKE_SUMMARY writes
 * smoke_summary.tsv with a GO / NO-GO verdict. No real data are used: genotype and
 * phenotype values are synthetic software fixtures generated inside the tasks.
 * Run it twice with -resume to check task caching (every task should be CACHED).
 */

process SMOKE_R {
    label 'process_single'
    errorStrategy 'ignore'

    output:
    path 'result.tsv'

    script:
    """
    rv=\$(Rscript -e 'cat(R.version.string)' 2>&1) || true
    if echo "\$rv" | grep -q "R version 4.3.1"; then st=PASS; else st=FAIL; fi
    printf "S01_container_R\\t%s\\t%s host=%s\\n" "\$st" "\$rv" "\$(hostname)" > result.tsv
    """
}

process SMOKE_TOOLS {
    label 'process_single'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/smoke", mode: params.publish_mode, pattern: 'tool_versions.tsv'

    output:
    path 'result.tsv', emit: result
    path 'tool_versions.tsv', emit: versions

    script:
    """
    xqtl_tool_versions.sh > tool_versions.tsv 2> tool_versions.err || true
    nmiss=\$(awk -F'\\t' 'NR>1 && \$2=="MISSING"' tool_versions.tsv | wc -l)
    ntot=\$(awk 'NR>1' tool_versions.tsv | wc -l)
    if [ "\$ntot" -gt 0 ] && [ "\$nmiss" -eq 0 ]; then st=PASS; else st=FAIL; fi
    printf "S02_embedded_tools\\t%s\\t%s tools checked, %s missing: %s\\n" "\$st" "\$ntot" "\$nmiss" \
        "\$(awk -F'\\t' 'NR>1 && \$2=="MISSING"{printf "%s ", \$1}' tool_versions.tsv)" > result.tsv
    """
}

process SMOKE_REFS {
    label 'process_single'
    errorStrategy 'ignore'

    output:
    path 'result.tsv'

    script:
    """
    if xqtl_refs_check.py --refs-dir '${params.refs_dir}' --quick > refs.log 2>&1; then st=PASS; else st=FAIL; fi
    printf "S03_embedded_references\\t%s\\t%s\\n" "\$st" "\$(tail -1 refs.log | tr '\\t' ' ')" > result.tsv
    """
}

process SMOKE_RESTRICTED {
    label 'process_single'
    errorStrategy 'ignore'

    input:
    val archive

    output:
    path 'result.tsv'

    script:
    if (!archive)
        """
        printf "S04_restricted_file\\tWARN\\tnot set (--restricted); request xqtl-restricted-<version>.tar.gz by email before site runs (docs/SETUP.md)\\n" > result.tsv
        """
    else
        """
        if xqtl_restricted.py verify '${archive}' --expected-version '${params.expected_restricted_version}' > r.log 2>&1; then st=PASS; else st=FAIL; fi
        printf "S04_restricted_file\\t%s\\t%s\\n" "\$st" "\$(head -1 r.log | tr '\\t' ' ')" > result.tsv
        """
}

process SMOKE_PLINK2 {
    label 'process_single'
    errorStrategy 'ignore'

    output:
    path 'result.tsv'

    script:
    """
    set +e
    plink2 --dummy 200 2000 acgt dosage-freq=0.5 --make-pgen --out dummy > log.txt 2>&1
    plink2 --pfile dummy --freq --hardy --out dummy >> log.txt 2>&1
    plink2 --pfile dummy --r2-unphased --ld-window-kb 1000 --out dummy >> log.txt 2>&1
    ver=\$(plink2 --version 2>&1 | head -1)
    if [ -s dummy.afreq ] && [ -s dummy.hardy ] && [ -s dummy.vcor ]; then st=PASS; else st=FAIL; fi
    printf "S05_plink2\\t%s\\t%s afreq_rows=%s vcor_rows=%s\\n" "\$st" "\$ver" \
        "\$(wc -l < dummy.afreq 2>/dev/null)" "\$(wc -l < dummy.vcor 2>/dev/null)" > result.tsv
    """
}

process SMOKE_RESOURCES {
    cpus 4
    memory '8 GB'
    time '30m'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/smoke", mode: params.publish_mode, pattern: 'sched_view.txt'
    beforeScript 'if [ -n "${LSB_JOBID:-}" ] && command -v bjobs >/dev/null 2>&1; then bjobs -l "$LSB_JOBID" > sched_view.txt 2>&1; elif [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null 2>&1; then scontrol show job "$SLURM_JOB_ID" > sched_view.txt 2>&1; elif [ -n "${PBS_JOBID:-}" ] && command -v qstat >/dev/null 2>&1; then qstat -f "$PBS_JOBID" > sched_view.txt 2>&1; elif [ -n "${JOB_ID:-}" ] && command -v qstat >/dev/null 2>&1; then qstat -j "$JOB_ID" > sched_view.txt 2>&1; else echo "no scheduler query available" > sched_view.txt; fi'

    output:
    path 'result.tsv', emit: result
    path 'sched_view.txt', emit: sched

    script:
    """
    set +e
    ncpu=\$(nproc)
    Rscript -e 'x <- numeric(4e9/8); x[] <- 1; cat(sum(x) > 0)' > alloc.txt 2>&1
    alloc=\$(tail -c 200 alloc.txt | tr '\\n' ' ')
    if grep -q TRUE alloc.txt; then st=PASS; else st=FAIL; fi
    if [ "\$ncpu" -gt 4 ]; then note="CPUs not confined (every tool is given task.cpus explicitly)"; else note="CPUs confined to request"; fi
    printf "S06_resources\\t%s\\tnproc=%s (requested 4); 4GB_alloc=%s; %s\\n" "\$st" "\$ncpu" "\$alloc" "\$note" > result.tsv
    """
}

process SMOKE_FANOUT {
    label 'process_single'
    errorStrategy 'ignore'

    input:
    val i

    output:
    path 'result.tsv'

    script:
    """
    printf "S07_fanout_task%s\\tPASS\\thost=%s\\n" "${i}" "\$(hostname)" > result.tsv
    """
}

process SMOKE_ARRAY {
    label 'process_single'
    array 5
    errorStrategy 'ignore'

    input:
    val i

    output:
    path 'result.tsv'

    script:
    """
    printf "S08_array_task%s\\tPASS\\thost=%s index=%s\\n" "${i}" "\$(hostname)" \
        "\${LSB_JOBINDEX:-\${SLURM_ARRAY_TASK_ID:-\${PBS_ARRAY_INDEX:-\${SGE_TASK_ID:-NA}}}}" > result.tsv
    """
}

process SMOKE_BIND {
    label 'process_single'
    errorStrategy 'ignore'

    input:
    val p

    output:
    path 'result.tsv'

    script:
    if (!p)
        """
        printf "S09_bind_read\\tSKIP\\tset --bind_test_path to a small readable file on the project filesystem\\n" > result.tsv
        """
    else
        """
        if [ -r "${p}" ]; then
            printf "S09_bind_read\\tPASS\\tread %s bytes from ${p} inside container\\n" "\$(wc -c < '${p}')" > result.tsv
        else
            printf "S09_bind_read\\tFAIL\\t${p} not visible inside container (add its filesystem to --bind_paths)\\n" > result.tsv
        fi
        """
}

// Attempt 1 requests 2 GB and holds params.smoke_mem_hog_gb (3 GB) for 90 s; where the
// scheduler enforces memory it is killed and attempt 2 (4 GB) passes. The hold matters
// on sites that enforce memory by periodic sampling (seadragon) rather than cgroups.
process SMOKE_OOM_RETRY {
    memory { 2.GB * task.attempt }
    cpus 1
    time '30m'
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }

    output:
    path 'result.tsv'

    script:
    def hog = params.smoke_mem_hog_gb
    """
    Rscript -e 'x <- numeric(${hog}e9/8); x[] <- 1; Sys.sleep(90); cat("allocated\\n")'
    if [ ${task.attempt} -eq 1 ]; then
        printf "S10_mem_enforcement_retry\\tWARN\\tattempt 1 held ${hog} GB under a 2 GB request without being killed: memory is not enforced by this scheduler, so set --max_memory conservatively\\n" > result.tsv
    else
        printf "S10_mem_enforcement_retry\\tPASS\\tattempt 1 killed over its 2 GB limit; attempt ${task.attempt} succeeded with ${task.memory}\\n" > result.tsv
    fi
    """
}

process SMOKE_TENSORQTL {
    label 'process_gpu_capable'
    cpus 4
    memory '8 GB'
    time '1h'
    errorStrategy 'ignore'

    output:
    path 'result.tsv', emit: result
    path 'cis_perm.tsv', optional: true, emit: perm

    script:
    def dev = params.use_gpu ? 'cuda' : 'cpu'
    """
    tensorqtl_smoke.py --device ${dev} --threads ${task.cpus} --test-id S11_tensorqtl_${dev} --out result.tsv \
      || printf "S11_tensorqtl_${dev}\\tFAIL\\tpython exited non-zero (see .command.err)\\n" > result.tsv
    """
}

// Storey q-values through the atlas file bridge (tensorQTL -> TSV -> R qvalue); no rpy2.
process SMOKE_QVALUE {
    label 'process_single'
    errorStrategy 'ignore'

    input:
    path cis_perm

    output:
    path 'result.tsv'

    script:
    """
    awk -F'\\t' 'NR==1{print "pval_beta"; next}{print \$2}' ${cis_perm} > in.tsv
    if Rscript \$(command -v compute_qvalues.R) in.tsv out.tsv > bridge.log 2>&1; then
        Rscript -e '
          a <- read.delim("${cis_perm}"); q <- read.delim("out.tsv")
          ok <- nrow(q) == nrow(a) && all(q\$qval >= 0 & q\$qval <= 1) && q\$qval[a\$phenotype_id == "feat0"] <= 0.05
          cat(sprintf("S12_storey_qvalue\\t%s\\tcompute_qvalues.R qvalue=%s R=%s n=%d planted_feature_qval=%.2e\\n",
                      if (isTRUE(ok)) "PASS" else "FAIL", as.character(packageVersion("qvalue")),
                      getRversion(), nrow(q), q\$qval[a\$phenotype_id == "feat0"]), file = "result.tsv")
        ' || printf "S12_storey_qvalue\\tFAIL\\tcheck script error (see .command.err)\\n" > result.tsv
    else
        printf "S12_storey_qvalue\\tFAIL\\tcompute_qvalues.R failed: %s\\n" "\$(tail -2 bridge.log | tr '\\n\\t' '  ')" > result.tsv
    fi
    """
}

process SMOKE_SUMMARY {
    label 'process_local'
    publishDir "${params.outdir}/smoke", mode: params.publish_mode

    input:
    path results, stageAs: 'r??/*'
    val expected

    output:
    path 'smoke_summary.tsv'

    script:
    def exp = expected.join(' ')
    """
    printf "test_id\\tstatus\\tdetail\\n" > body.tsv
    cat r*/result.tsv | sort >> body.tsv
    for t in ${exp}; do
        grep -q "^\$t" body.tsv || printf "%s\\tFAIL\\tno result (task failed; see the trace and .command.err)\\n" "\$t" >> body.tsv
    done
    nfail=\$(awk -F'\\t' 'NR>1 && \$2=="FAIL"' body.tsv | wc -l)
    if [ "\$nfail" -eq 0 ]; then verdict=GO; else verdict="NO-GO (\$nfail FAIL)"; fi
    { cat body.tsv
      printf "S99_context\\tINFO\\tprofile=${workflow.profile} engine=${workflow.containerEngine ?: 'none'} nextflow=${nextflow.version} image=${params.sif ?: params.image}\\n"
      printf "VERDICT\\t%s\\t\\n" "\$verdict"; } > smoke_summary.tsv
    column -t -s \$'\\t' smoke_summary.tsv || cat smoke_summary.tsv
    """
}

workflow SMOKE {
    main:
    SMOKE_R()
    SMOKE_TOOLS()
    SMOKE_REFS()
    SMOKE_RESTRICTED(channel.value(params.restricted ?: ''))
    SMOKE_PLINK2()
    SMOKE_RESOURCES()
    SMOKE_FANOUT(channel.of(1..5))
    SMOKE_ARRAY(channel.of(1..5))
    SMOKE_BIND(channel.value(params.bind_test_path ?: ''))
    SMOKE_OOM_RETRY()
    SMOKE_TENSORQTL()
    SMOKE_QVALUE(SMOKE_TENSORQTL.out.perm)

    expected = ['S01_container_R', 'S02_embedded_tools', 'S03_embedded_references', 'S04_restricted_file',
                'S05_plink2', 'S06_resources', 'S07_fanout_task', 'S08_array_task', 'S09_bind_read',
                'S10_mem_enforcement_retry', 'S11_tensorqtl', 'S12_storey_qvalue']

    all_results = SMOKE_R.out
        .mix(SMOKE_TOOLS.out.result, SMOKE_REFS.out, SMOKE_RESTRICTED.out, SMOKE_PLINK2.out,
             SMOKE_RESOURCES.out.result, SMOKE_FANOUT.out, SMOKE_ARRAY.out, SMOKE_BIND.out,
             SMOKE_OOM_RETRY.out, SMOKE_TENSORQTL.out.result, SMOKE_QVALUE.out)
        .collect()

    SMOKE_SUMMARY(all_results, expected)

    emit:
    summary = SMOKE_SUMMARY.out
}
