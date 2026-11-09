/*
 * --mode imputation_prep : pre-imputation genotypes -> per-chromosome GRCh38 VCFs for the TOPMed
 * Imputation Server (atlas 01_genotype_imputation: liftover + filter_vcf_for_imputation).
 * Scripts: bin/xqtl_impprep.sh (convert, filter), bin/xqtl_impprep.py (sex check, report).
 */

process IP_GENOME {
    tag "${fasta.name}"
    label 'process_low'

    input:
    path fasta
    path fasta_idx, stageAs: 'idx/*'
    val build

    output:
    tuple path('genome'), val(fasta.name), emit: genome

    script:
    def dict = fasta.name.replaceAll(/\.gz$/, '').replaceAll(/\.(fa|fasta|fna)$/, '') + '.dict'
    """
    mkdir genome
    ln -s ../${fasta} genome/${fasta.name}
    for f in idx/*; do [ -e "\$f" ] && ln -s ../\$f genome/\$(basename \$f); done
    [ -e genome/${fasta.name}.fai ] || samtools faidx genome/${fasta.name}
    if [ "${build}" = GRCh37 ]; then
        picard CreateSequenceDictionary R=genome/${fasta.name} O=genome/${dict} > dict.log 2>&1
    fi
    """
}

process IP_CONVERT {
    tag "${kind} ${prefix}"
    label 'process_medium'

    input:
    tuple val(kind), path(raw, stageAs: 'raw/*'), val(prefix), val(build)

    output:
    tuple path('converted.vcf.gz'), path('converted.vcf.gz.csi'), emit: vcf
    path 'chroms.tsv', emit: chroms
    path 'sample_missing.smiss', emit: smiss
    path 'sexcheck.tsv', emit: sexcheck

    script:
    """
    xqtl_impprep.sh convert ${kind} raw/${prefix} ${build} ${task.cpus}
    xqtl_impprep.py sexcheck --gt chrX_gt.tsv --samples chrX_samples.txt --out sexcheck.tsv
    """
}

process IP_LIFTOVER {
    tag 'GRCh37 -> GRCh38'
    label 'process_medium'

    input:
    tuple path(vcf), path(idx)
    tuple path(genome), val(fa)
    path chain

    output:
    tuple path('lifted.vcf.gz'), path('lifted.vcf.gz.tbi'), emit: vcf
    path 'chroms_lifted.tsv', emit: chroms
    path 'liftover.tsv', emit: counts

    script:
    def mem = Math.max(2, (task.memory.toGiga() * 0.8) as int)
    """
    picard -Xmx${mem}g LiftoverVcf I=${vcf} O=lifted.vcf.gz CHAIN=${chain} REJECT=rejected.vcf.gz \\
        R=${genome}/${fa} WARN_ON_MISSING_CONTIG=true RECOVER_SWAPPED_REF_ALT=true CREATE_INDEX=false \\
        TMP_DIR=. > liftover.log 2>&1 || { tail -20 liftover.log >&2; exit 1; }
    # htslib index (htsjdk .tbi files lack the record counts used by bcftools index -n/-s)
    tabix -f -p vcf lifted.vcf.gz
    xqtl_impprep.py liftover --rejected rejected.vcf.gz --lifted \$(bcftools index -n lifted.vcf.gz) --out liftover.tsv
    bcftools index -s lifted.vcf.gz | awk 'BEGIN{OFS="\\t"; print "chrom", "records"} {print \$1, \$3}' > chroms_lifted.tsv
    """
}

process IP_FILTER {
    tag "${chrom}"
    label 'process_low'

    input:
    tuple val(chrom), path(vcf), path(idx)
    tuple path(genome), val(fa)
    val ip

    output:
    tuple val(chrom), path("${chrom}.vcf.gz"), path("${chrom}.vcf.gz.tbi"), emit: vcf
    path "${chrom}.counts.tsv", emit: counts

    script:
    """
    xqtl_impprep.sh filter ${vcf} ${chrom} ${genome}/${fa} ${ip.geno_missing_max} ${ip.maf_min} ${ip.hwe_min} ${task.cpus}
    """
}

process IP_REPORT {
    label 'process_single'
    publishDir path: "${params.outdir}/imputation_prep", mode: params.publish_mode

    input:
    path counts, stageAs: 'counts/*'
    path smiss
    path sexcheck
    path liftover, stageAs: 'lift/*'
    path metadata, stageAs: 'meta/*'
    val cfg

    output:
    path 'imputation_prep_report.tsv', emit: report
    path 'imputation_prep_samples.tsv', emit: samples
    path 'TOPMED_CHECKLIST.md', emit: checklist
    path 'liftover_summary.tsv', optional: true, emit: liftover

    script:
    def lift = liftover ? "--liftover ${liftover}" : ''
    def meta = metadata ? "--metadata ${metadata}" : ''
    """
    cat > cfg.json <<'JSON'
${groovy.json.JsonOutput.toJson(cfg)}
JSON
    xqtl_impprep.py report --counts counts/*.counts.tsv --smiss ${smiss} --sexcheck ${sexcheck} ${lift} ${meta} \\
        --config cfg.json --out-dir .
    ${liftover ? "cp ${liftover} liftover_summary.tsv" : ''}
    """
}

process IP_PUBLISH {
    label 'process_local'
    publishDir path: "${params.outdir}/imputation_prep", mode: params.publish_mode

    input:
    path vcfs, stageAs: 'vcf/*'
    path report

    output:
    path 'chr*.vcf.gz*', emit: vcfs

    script:
    """
    for f in vcf/*; do ln -s \$f \$(basename \$f); done
    """
}
