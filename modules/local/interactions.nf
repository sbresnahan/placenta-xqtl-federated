/*
 * Interaction processes (M5): signal test sets, G x C at lead signals (atlas 49), G x E exposures,
 * tier-1 / tier-2 scans sharded by chromosome, merges and the SNP x covariate sensitivity
 * (atlas 07_gxe_mapping 50/52/54/55). Outputs are published to ${outdir}/site/<ANC>/ with the
 * <ANC>_ prefix that EXPORT_STRATUM and --mode export read; the per-individual exposure table and
 * the extracted MDACC signal list stay in the work directory.
 */

process RESTRICTED_SIGNALS {
    label 'process_local'

    input:
    path archive
    val expected

    output:
    path 'mdacc_signals.tsv.gz', emit: signals

    script:
    def ev = expected ? "--expected-version ${expected}" : ''
    """
    xqtl_restricted.py extract ${archive} restricted ${ev}
    mv restricted/signals/mdacc_signals.tsv.gz mdacc_signals.tsv.gz
    rm -rf restricted
    """
}

process GX_SIGNALS {
    tag "${anc} ${mod}"
    label 'process_low'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode

    input:
    tuple val(anc), val(mod), path(bed), path(cis), path(pgen), path(pvar), path(psam)
    path signals
    val m

    output:
    tuple val(anc), val(mod), path("${anc}_${mod}.gx_pairs.tsv"), path("${anc}_${mod}.gx_phenotypes.tsv"), emit: sets

    script:
    """
    xqtl_signals.py --anc ${anc} --modality ${mod} --cis ${cis} --bed ${bed} --pvar ${pvar} \\
        --signals ${signals} --fdr ${m.fdr} --out ${anc}_${mod}
    """
}

process GXC {
    tag "${anc} ${mod}"
    label 'process_medium'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode

    input:
    tuple val(anc), val(mod), path(bed), path(pairs), path(cov), path(samples), path(pgen), path(pvar), path(psam)
    path cell_props
    val g

    output:
    tuple val(anc), val(mod), path("${anc}_${mod}.gxc.parquet"), path("${anc}_${mod}.gxc_annotation.tsv.gz"), emit: gxc

    script:
    def ex = "--exclude ${(g.exclude_cell_types ?: []).join(' ')}"
    """
    xqtl_gxc.py --anc ${anc} --modality ${mod} --pairs ${pairs} --bed ${bed} --cov ${cov} --pfile ${pgen.baseName} \\
        --samples ${samples} --cell-props ${cell_props} ${ex} --min-n ${g.min_n} \\
        --maternal-threshold ${g.maternal_threshold} --fdr ${g.fdr} --workers ${task.cpus} --out ${anc}_${mod}
    """
}

process GXE_EXPOSURES {
    label 'process_single'
    publishDir path: "${params.outdir}/site", mode: params.publish_mode,
        saveAs: { fn -> fn == 'gxe_exposures.tsv' ? null : "${fn.tokenize('.')[0]}/${fn}" }

    input:
    path metadata
    tuple val(ancs), path(samples, stageAs: 'smp/*')
    val g

    output:
    path 'gxe_exposures.tsv', emit: table
    path '*.gxe_exposures.tsv', emit: summary

    script:
    def js = groovy.json.JsonOutput.toJson(g)
    def smp = [ancs, samples].transpose().collect { a, f -> "${a}=${f}" }.join(' ')
    """
    cat > gxe.json <<'JSON'
    ${js}
    JSON
    xqtl_gxe.py exposures --metadata ${metadata} --config gxe.json --samples ${smp}
    """
}

process GXE_SCAN {
    tag "${anc} ${mod} ${exp} tier ${tier} ${chrom}"
    label 'process_gpu_capable'

    input:
    tuple val(anc), val(mod), val(exp), val(tier), val(chrom), path(bed), path(cov), path(pgen), path(pvar), path(psam), path(summary), path(phenos)
    path exposures
    val g

    output:
    tuple val(anc), val(mod), val(exp), val(tier), path("${anc}_${mod}_${exp}.t${tier}.${chrom}.parquet"),
        path("${anc}_${mod}_${exp}.t${tier}.${chrom}.untested.tsv"), path("${anc}_${mod}_${exp}.t${tier}.${chrom}.status.tsv"), emit: shard

    script:
    def pre = "${anc}_${mod}_${exp}.t${tier}.${chrom}"
    def ph = tier == 2 ? "--phenotypes ${phenos} --tier2-variants ${g.tier2}" : ''
    def ds = g.dosages ? '--dosages' : ''
    """
    xqtl_gxe.py scan --tier ${tier} --anc ${anc} --modality ${mod} --exposure ${exp} --chrom ${chrom} \\
        --bed ${bed} --cov ${cov} --pfile ${pgen.baseName} --exposures ${exposures} --summary ${summary} ${ph} \\
        --maf-interaction ${g.maf_interaction} --min-n ${g.min_n} --blocks ${g.tier1_permutation_blocks.join(' ')} \\
        --seed ${g.seed} ${ds} --threads ${task.cpus} --out ${pre}
    """
}

process GXE_MERGE {
    tag "${anc} ${mod} ${exp} tier ${tier}"
    label 'process_low'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode

    input:
    tuple val(anc), val(mod), val(exp), val(tier), path(shards, stageAs: 'shards/*'), path(untested, stageAs: 'shards/*'), path(status, stageAs: 'shards/*'), path(phenos)
    val g

    output:
    tuple val(anc), val(mod), val(exp), val(tier), path("${anc}_${mod}_${exp}.gxe_tier${tier}.parquet"), emit: merged
    tuple val(anc), path("${anc}_${mod}_${exp}.gxe_tier*"), emit: files

    script:
    def ph = tier == 2 ? "--phenotypes ${phenos}" : ''
    """
    xqtl_gxe.py merge --tier ${tier} --anc ${anc} --modality ${mod} --exposure ${exp} --shards ${shards} ${ph} \\
        --fdr ${g.fdr} --out ${anc}_${mod}_${exp}
    """
}

process GXE_SENSITIVITY {
    tag "${anc} ${mod} ${exp}"
    label 'process_low'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode

    input:
    tuple val(anc), val(mod), val(exp), path(tier1), path(bed), path(cov), path(samples), path(pgen), path(pvar), path(psam)
    path exposures
    path metadata
    val g

    output:
    tuple val(anc), path("${anc}_${mod}_${exp}.gxe_sensitivity.tsv"), path("${anc}_${mod}_${exp}.gxe_sensitivity.status.tsv"), emit: files

    script:
    def js = groovy.json.JsonOutput.toJson(g)
    """
    cat > gxe.json <<'JSON'
    ${js}
    JSON
    xqtl_gxe.py sensitivity --anc ${anc} --modality ${mod} --exposure ${exp} --tier1 ${tier1} --bed ${bed} \\
        --cov ${cov} --pfile ${pgen.baseName} --exposures ${exposures} --metadata ${metadata} --samples ${samples} \\
        --config gxe.json --fdr ${g.fdr} --out ${anc}_${mod}_${exp}
    """
}
