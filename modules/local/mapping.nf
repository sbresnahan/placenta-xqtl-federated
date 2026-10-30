/*
 * Mapping processes: normalization within stratum, HCP factor grid, k selection, tensorQTL
 * permutation / conditionally independent / nominal passes and the combined-modality arm.
 * Staging files are published to ${outdir}/site/<ANC>/ with the names --mode export expects.
 */

process NORMALIZE {
    tag "${anc} ${mod}"
    label 'process_low'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode,
        pattern: '*.{features.tsv,norm_summary.tsv,outliers.tsv}'

    input:
    tuple val(anc), val(mod), path(bed), path(counts), path(exclude), path(samples)

    output:
    tuple val(anc), val(mod), path("${anc}_${mod}.bed.gz"), path("${anc}_${mod}.phenotype_groups.txt", arity: '0..*'), emit: bed
    tuple val(anc), val(mod), path("${anc}_${mod}.features.tsv"), path("${anc}_${mod}.norm_summary.tsv"), emit: staging
    tuple val(anc), val(mod), path("${anc}_${mod}.outliers.tsv"), emit: outliers

    script:
    def c = counts ? "--counts ${counts}" : ''
    def x = exclude ? "--exclude ${exclude}" : ''
    """
    xqtl_normalize.R --bed ${bed} ${c} ${x} --modality ${mod} --samples ${samples} --out ${anc}_${mod}
    """
}

process COMBINE_MODALITIES {
    tag "${anc}"
    label 'process_low'

    input:
    tuple val(anc), val(mods), path(beds, stageAs: 'in/*')

    output:
    tuple val(anc), val('combined'), path("${anc}_combined.bed.gz"), path("${anc}_combined.phenotype_groups.txt"), emit: bed

    script:
    def args = [mods, beds].transpose().collect { m, b -> "${m}=${b}" }.join(' ')
    """
    xqtl_combine_modalities.py --out ${anc}_combined ${args}
    """
}

// HCP factors at one k, the full covariate table, and the eGene count used to pick k*.
process HCP_EVAL {
    tag "${anc} ${arm} k=${k}"
    label 'process_gpu_capable'

    input:
    tuple val(anc), val(arm), val(k), path(bed), path(groups), path(samples), path(pcs), path(pgen), path(pvar), path(psam)
    path metadata
    path qc_metrics
    path cell_props
    val m

    output:
    tuple val(anc), val(arm), val(k), path("${anc}_${arm}.k${k}.eval.tsv"), path("${anc}_${arm}.cov_k${k}.tsv"), emit: eval

    script:
    def g = groups ? "--groups ${groups}" : ''
    def ct = cell_props ? "--cell-props ${cell_props}" : ''
    def ex = m.exclude_cell_types ? "--exclude ${m.exclude_cell_types.join(' ')}" : ''
    """
    xqtl_hcp.R --bed ${bed} --qc-metrics ${qc_metrics} --samples ${samples} --k ${k} \\
        --max-phenotypes ${m.hcp_max_phenotypes} --seed ${m.seed} --out ${anc}_${arm}.hcp_k${k}.tsv
    xqtl_covariates.py --bed ${bed} --samples ${samples} --metadata ${metadata} --geno-pcs ${pcs} ${ct} \\
        --hcp ${anc}_${arm}.hcp_k${k}.tsv ${ex} --out ${anc}_${arm}.cov_k${k}.tsv
    xqtl_tensorqtl.py eval --bed ${bed} --geno ${pgen.baseName} --cov ${anc}_${arm}.cov_k${k}.tsv ${g} \\
        --nperm ${m.nperm} --window ${m.cis_window} --maf ${m.maf} --fdr ${m.fdr} --seed ${m.seed} \\
        --eval-chrom ${m.eval_chrom} --eval-min ${m.eval_min_phenotypes} --k ${k} --threads ${task.cpus} \\
        --out ${anc}_${arm}.k${k}
    """
}

process SELECT_K {
    tag "${anc} ${mod}"
    label 'process_local'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode, pattern: '*.settings.json'

    input:
    tuple val(anc), val(mod), val(arms), val(grouped), path(norm_summary), path(evals)
    val m

    output:
    tuple val(anc), val(mod), path('k.txt'), emit: k
    tuple val(anc), path('*.settings.json', arity: '1..*'), emit: settings

    script:
    def ns = norm_summary ? "--norm-summary ${norm_summary}" : ''
    def armargs = arms.collect { a -> "--arm ${a}" }.join(' ')
    """
    xqtl_select_k.py --anc ${anc} --modality ${mod} ${armargs} ${grouped ? '--grouped' : ''} ${ns} \\
        --nperm ${m.nperm} --window ${m.cis_window} --maf ${m.maf} --fdr ${m.fdr} --seed ${m.seed} \\
        --out-k k.txt ${evals}
    """
}

process CIS_PERMUTATION {
    tag "${anc} ${arm}"
    label 'process_gpu_capable'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode

    input:
    tuple val(anc), val(arm), path(bed), path(groups), path(cov), path(pgen), path(pvar), path(psam)
    val m

    output:
    tuple val(anc), val(arm), path("${anc}_${arm}.cis.parquet"), path("${anc}_${arm}.cis_summary.tsv"), emit: cis
    tuple val(anc), val(arm), path("${anc}_${arm}.covariates.tsv"), emit: cov

    script:
    def g = groups ? "--groups ${groups}" : ''
    """
    [ "${cov}" = "${anc}_${arm}.covariates.tsv" ] || cp ${cov} ${anc}_${arm}.covariates.tsv
    xqtl_tensorqtl.py cis --bed ${bed} --geno ${pgen.baseName} --cov ${anc}_${arm}.covariates.tsv ${g} \\
        --nperm ${m.nperm} --window ${m.cis_window} --maf ${m.maf} --fdr ${m.fdr} --seed ${m.seed} \\
        --threads ${task.cpus} --out ${anc}_${arm}
    """
}

process CIS_INDEPENDENT {
    tag "${anc} ${arm} chunk ${chunk}"
    label 'process_gpu_capable'

    input:
    tuple val(anc), val(arm), val(chunk), path(bed), path(groups), path(cov), path(pgen), path(pvar), path(psam), path(cis)
    val m

    output:
    tuple val(anc), val(arm), path("chunk_*.parquet"), emit: chunk

    script:
    def g = groups ? "--groups ${groups}" : ''
    def name = String.format('chunk_%05d', chunk as int)
    """
    xqtl_tensorqtl.py independent --bed ${bed} --geno ${pgen.baseName} --cov ${cov} ${g} --cis-table ${cis} \\
        --chunk-index ${chunk} --chunk-size ${m.independent_chunk_size} \\
        --nperm ${m.nperm} --window ${m.cis_window} --maf ${m.maf} --fdr ${m.fdr} --seed ${m.seed} \\
        --threads ${task.cpus} --out ${name}
    """
}

process MERGE_INDEPENDENT {
    tag "${anc} ${arm}"
    label 'process_local'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode

    input:
    tuple val(anc), val(arm), path(chunks, stageAs: 'chunks/*')

    output:
    tuple val(anc), val(arm), path("${anc}_${arm}.independent.parquet"), emit: independent

    script:
    """
    xqtl_merge_parquet.py ${anc}_${arm}.independent.parquet \$(ls chunks/*.parquet 2>/dev/null | sort)
    """
}

process CIS_NOMINAL {
    tag "${anc} ${mod}"
    label 'process_gpu_capable'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode

    input:
    tuple val(anc), val(mod), path(bed), path(cov), path(pgen), path(pvar), path(psam)
    val m

    output:
    tuple val(anc), val(mod), path("${anc}_${mod}.nominal"), emit: nominal

    script:
    """
    xqtl_tensorqtl.py nominal --bed ${bed} --geno ${pgen.baseName} --cov ${cov} \\
        --window ${m.cis_window} --maf ${m.maf} --threads ${task.cpus} --out ${anc}_${mod}
    """
}
