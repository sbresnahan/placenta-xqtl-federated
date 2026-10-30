/*
 * Export processes: assemble files A-F and I per stratum, then I_strata, MANIFEST (J) and the
 * data dictionary for the cohort, followed by the mandatory leak scan. Only ${outdir}/export/
 * is meant to leave the site.
 */

process VERIFY_RESTRICTED {
    label 'process_local'
    publishDir "${params.outdir}/site", mode: params.publish_mode

    input:
    path archive
    val expected

    output:
    path 'restricted_version.txt', emit: version

    script:
    def ev = expected ? "--expected-version ${expected}" : ''
    """
    xqtl_restricted.py verify ${archive} ${ev} | tee verify.log
    sed -n 's/^restricted file OK: version \\([^;]*\\);.*/\\1/p' verify.log > restricted_version.txt
    [ -s restricted_version.txt ]
    """
}

process EXPORT_STRATUM {
    tag "${anc}"
    label 'process_low'

    input:
    tuple val(anc), path(staged, stageAs: 'stage/*')
    path qc_metrics
    path cell_props
    val ex

    output:
    path "${anc}", type: 'dir', emit: dir

    script:
    def ct = cell_props ? "--cell-props ${cell_props}" : ''
    def rd = ex.round_digits != null ? "--round-digits ${ex.round_digits}" : ''
    """
    xqtl_export.py stratum --anc ${anc} --cohort ${ex.cohort_id} --in stage --qc-metrics ${qc_metrics} ${ct} \\
        --format ${ex.format} --min-mac ${ex.min_mac_export} ${rd} --out ${anc}
    """
}

process EXPORT_FINALIZE {
    label 'process_low'
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'export'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_mode, pattern: 'leak_report.tsv'

    input:
    path strata_dirs, stageAs: 'strata/*'
    path strata
    path ids, stageAs: 'ids/*'
    val ex

    output:
    path 'export', type: 'dir', emit: export
    path 'leak_report.tsv', emit: leak_report

    script:
    """
    xqtl_export.py finalize --cohort ${ex.cohort_id} --strata ${strata} --ids ${ids} \\
        --image-sha256 ${ex.image_sha256} --restricted-version ${ex.restricted_version} \\
        --pipeline-version ${ex.pipeline_version} --leak-report leak_report.tsv --out export ${strata_dirs}
    """
}
