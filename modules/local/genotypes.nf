/*
 * Genotype processes: post-imputation QC, merge, ancestry assignment, per-stratum QC and LD.
 * Everything written here stays at the site; only the export step (modules/local/export.nf)
 * produces shareable files.
 */

process GENO_QC {
    tag "${chrom}"
    label 'process_low'

    input:
    tuple val(chrom), path(vcf)
    path metadata
    val rsq_min

    output:
    path "${chrom}.*", emit: files

    script:
    """
    awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="genotype_id") c=i; next} c && \$c!="" && \$c!="NA"{print \$c}' \\
        ${metadata} | sort -u > genotype_ids.txt
    xqtl_geno_qc.sh --vcf ${vcf} --samples genotype_ids.txt --chrom ${chrom} --rsq-min ${rsq_min} \\
        --out ${chrom} --threads ${task.cpus}
    """
}

process GENO_MERGE {
    label 'process_low'
    publishDir "${params.outdir}/site/genotype_qc", mode: params.publish_mode, pattern: 'cohort.{qc_counts.tsv,vcf_samples.txt}'

    input:
    path files
    val chroms

    output:
    tuple path('cohort.pgen'), path('cohort.pvar'), path('cohort.psam'), emit: pfile
    path 'cohort.info.tsv', emit: info
    path 'cohort.qc_counts.tsv', emit: counts
    path 'cohort.vcf_samples.txt', emit: vcf_samples

    script:
    """
    xqtl_geno_merge.sh --out cohort --threads ${task.cpus} ${chroms.join(' ')}
    """
}

process ANCESTRY {
    label 'process_low'
    publishDir "${params.outdir}/site/ancestry", mode: params.publish_mode

    input:
    tuple path(pgen), path(pvar), path(psam)
    path panel, stageAs: 'panel/*'
    path superpop
    path metadata
    val anc_cfg

    output:
    path 'cohort_strata.tsv', emit: strata
    path 'cohort_ancestry.tsv', emit: ancestry
    path 'stratum_*.samples.tsv', emit: samples, optional: true
    path 'cohort_pca.png', emit: png
    path 'cohort_pcair_pcs*.tsv', emit: pcs

    script:
    """
    xqtl_ancestry.sh --pfile cohort --panel panel/1kGP_hg38_pca_panel --superpop ${superpop} \\
        --metadata ${metadata} --min-n ${anc_cfg.min_n} --prob-min ${anc_cfg.assign_prob_min} \\
        --groups ${anc_cfg.groups.join(',')} --out cohort --threads ${task.cpus}
    """
}

process STRATUM_GENO {
    tag "${anc}"
    label 'process_low'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode,
        pattern: '*.{variant_qc.tsv,geno_qc_counts.tsv,samples.tsv,geno_outliers.tsv}'

    input:
    tuple val(anc), path(samples), path(pgen), path(pvar), path(psam), path(info)
    val qc

    output:
    tuple val(anc), path("${anc}.pgen"), path("${anc}.pvar"), path("${anc}.psam"), emit: geno
    tuple val(anc), path("${anc}.geno_pcs.tsv"), emit: pcs
    tuple val(anc), path("${anc}.samples.tsv"), emit: samples
    tuple val(anc), path("${anc}.variant_qc.tsv"), emit: variant_qc
    tuple val(anc), path("${anc}.variant_qc.tsv"), path("${anc}.geno_qc_counts.tsv"), path("${anc}.samples.tsv"), emit: staging
    path "${anc}.geno_outliers.tsv", emit: outliers

    script:
    """
    xqtl_stratum_geno.sh --pfile ${pgen.baseName} --info ${info} --samples ${samples} --anc ${anc} \\
        --maf ${qc.maf_min} --mac ${qc.mac_min} --hwe ${qc.hwe_p_min} --geno ${qc.geno_missing_max} \\
        --mind ${qc.sample_missing_max} --threads ${task.cpus}
    """
}

process LD_MATRICES {
    tag "${anc}"
    label 'process_low'
    publishDir path: { "${params.outdir}/site/${anc}" }, mode: params.publish_mode

    input:
    tuple val(anc), path(pgen), path(pvar), path(psam), path(variant_qc)
    val ld

    output:
    tuple val(anc), path("${anc}_ld"), emit: ld

    script:
    """
    xqtl_ld.py --geno ${anc} --variant-qc ${variant_qc} --out ${anc}_ld \\
        --maf-min ${ld.maf_min} --rsq-min ${ld.rsq_min} --window-mb ${ld.window_mb} --step-mb ${ld.step_mb}
    """
}
