/*
 * Phenotyping processes (Pantry-based): per-sample alignment and
 * quantification, cohort-level assembly of the nine phenotype BEDs, Picard QC metrics and MuSiC
 * cell-type proportions. Tool settings are listed in docs/PHENOTYPING.md.
 * Per-sample intermediates stay in the work directory; finished tables are published by
 * PHENO_FINALIZE to ${outdir}/phenotypes/.
 */

def asList(x) {
    return x instanceof List ? x : [x]
}

// ---------------------------------------------------------------- per sample ----------------

process STAR_ALIGN {
    tag "${sample}"
    label 'process_high'

    input:
    tuple val(sample), path(reads, stageAs: 'reads/*')
    path star

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam
    tuple val(sample), path("${sample}.Log.final.out"), emit: log

    script:
    def fq = asList(reads)
    def rfc = fq[0].name.endsWith('.gz') ? '--readFilesCommand gunzip -c' : ''
    """
    STAR --runMode alignReads --runThreadN ${task.cpus} --genomeDir ${star} \\
        --readFilesIn ${fq.join(' ')} ${rfc} \\
        --twopassMode Basic --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \\
        --outFileNamePrefix ${sample}.
    mv ${sample}.Aligned.sortedByCoord.out.bam ${sample}.bam
    samtools index ${sample}.bam
    rm -rf ${sample}._STARgenome ${sample}._STARpass1 ${sample}._STARtmp
    """
}

// expression: annotation transcripts, 20 bootstraps (needed by the QU correction)
process SALMON_EXPRESSION {
    tag "${sample}"
    label 'process_medium'

    input:
    tuple val(sample), path(reads, stageAs: 'reads/*')
    path index

    output:
    tuple val(sample), path("${sample}"), emit: quant

    script:
    def fq = asList(reads)
    def rd = fq.size() == 2 ? "-1 ${fq[0]} -2 ${fq[1]}" : "-r ${fq[0]}"
    """
    salmon quant -i ${index} -l A --validateMappings --numBootstraps 20 ${rd} -p ${task.cpus} -o ${sample}
    test -s ${sample}/quant.sf
    """
}

// alternative TSS / polyA: txrevise upstream and downstream event sets (no bootstraps)
process SALMON_TXREVISE {
    tag "${sample}"
    label 'process_medium'

    input:
    tuple val(sample), path(reads, stageAs: 'reads/*')
    path indices, stageAs: 'idx/*'

    output:
    tuple val(sample), path("${sample}.txrevise"), emit: quant

    script:
    def fq = asList(reads)
    def rd = fq.size() == 2 ? "-1 ${fq[0]} -2 ${fq[1]}" : "-r ${fq[0]}"
    """
    for idx in idx/*; do
        set_name=\$(basename "\$idx")
        salmon quant -i "\$idx" -l A --validateMappings ${rd} -p ${task.cpus} -o "${sample}.txrevise/\$set_name"
        test -s "${sample}.txrevise/\$set_name/quant.sf"
    done
    """
}

process REGTOOLS_JUNCTIONS {
    tag "${sample}"
    label 'process_single'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path "${sample}.junc", emit: junc

    script:
    """
    regtools junctions extract -a 8 -m 50 -M 500000 -s XS ${bam} -o ${sample}.junc
    """
}

// exon / intron read counts for the stability modality
process FEATURECOUNTS {
    tag "${sample}"
    label 'process_low'

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(exonic), path(intronic)
    val rnaseq

    output:
    path "${sample}.{exonic,intronic}.counts.txt", emit: counts

    script:
    def s = [reverse: 2, forward: 1, unstranded: 0][rnaseq.strandedness]
    def p = rnaseq.read_layout == 'paired' ? '-p' : ''
    """
    featureCounts ${bam} ${p} -a ${exonic} -t exon -s ${s} --fracOverlap 1 -T ${task.cpus} -o ${sample}.exonic.counts.txt
    featureCounts ${bam} ${p} -a ${intronic} -t intron -s ${s} --fracOverlap 0 -T ${task.cpus} -o ${sample}.intronic.counts.txt
    """
}

process EDITING_LEVELS {
    tag "${sample}"
    label 'process_single'

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(fa), path(fai), path(dict)
    tuple path(sites), path(site_genes)

    output:
    path "${sample}.rnaeditlevel.tsv.gz", emit: levels

    script:
    """
    query_editing_level.py --edit_sites ${sites} --ref ${fa} --bam ${bam} --output ${sample}.rnaeditlevel.tsv.gz
    """
}

process IRFINDER {
    tag "${sample}"
    label 'process_low'

    input:
    tuple val(sample), path(bam), path(bai)
    path irf_ref, stageAs: 'REF/IRFinder'

    output:
    path "${sample}.IRFinder-IR-*dir.txt", emit: ir

    script:
    """
    IRFinder BAM -r REF -d irf ${bam}
    for f in irf/IRFinder-IR-dir.txt irf/IRFinder-IR-nondir.txt; do
        [ -s "\$f" ] && cp "\$f" "${sample}.\$(basename "\$f")"
    done
    ls ${sample}.IRFinder-IR-*dir.txt > /dev/null
    """
}

process PICARD_QC {
    tag "${sample}"
    label 'process_low'

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(fa), path(fai), path(dict)
    path refflat
    val rnaseq

    output:
    path "${sample}.qcmetrics.tsv", emit: row

    script:
    def p = rnaseq.read_layout == 'paired' ? '--paired' : ''
    def mem = Math.max(1, (task.memory.toGiga() * 0.8) as int)
    """
    xqtl_picard_qc.py collect --bam ${bam} --sample ${sample} --refflat ${refflat} --fasta ${fa} \\
        --outdir picard_raw --output ${sample}.qcmetrics.tsv ${p} --java-mem ${mem}g
    """
}

// ---------------------------------------------------------------- cohort level ---------------

process AGG_EXPRESSION {
    label 'process_medium'

    input:
    path quants, stageAs: 'salmon/*'
    path samples
    path gtf

    output:
    path 'expression.bed', emit: expression
    path 'isoform_expression.bed', emit: isoform_expression
    path 'isoform_usage.bed', emit: isoform_usage
    path '{expression,isoform_expression}.counts.bed', emit: counts

    script:
    """
    xqtl_assemble_bed.py expression --samples ${samples} --input-dir salmon --ref-anno ${gtf} \\
        --output-expression expression.bed --skip-feature-filter
    xqtl_qu_correct_salmon.R --samples ${samples} --salmon-dir salmon --ref-anno ${gtf} --out-dir qu --chunk-size 10
    xqtl_assemble_bed.py expression --samples ${samples} --input-dir qu --ref-anno ${gtf} --units tpm_from_counts \\
        --output-isoforms isoform_usage.bed --output-isoform-expr isoform_expression.bed --skip-feature-filter
    xqtl_tximport_counts.R --samples ${samples} --salmon-dir salmon --qu-dir qu --ref-anno ${gtf} --out-dir .
    mv expression_counts.bed expression.counts.bed
    mv isoform_expression_counts.bed isoform_expression.counts.bed
    """
}

process AGG_ALT_ENDS {
    label 'process_low'

    input:
    path quants
    path samples
    path gtf

    output:
    path 'alt_tss.bed', emit: alt_tss
    path 'alt_polya.bed', emit: alt_polya

    script:
    """
    for s in \$(cat ${samples}); do
        for set_name in grp_1.upstream grp_2.upstream grp_1.downstream grp_2.downstream; do
            mkdir -p "\$set_name"
            ln -s "../\$s.txrevise/txrevise.\$set_name" "\$set_name/\$s"
        done
    done
    xqtl_assemble_bed.py alt-tss-polya --samples ${samples} --group1-dir grp_1.upstream --group2-dir grp_2.upstream \\
        --ref-anno ${gtf} --output alt_tss.bed --skip-feature-filter
    xqtl_assemble_bed.py alt-tss-polya --samples ${samples} --group1-dir grp_1.downstream --group2-dir grp_2.downstream \\
        --ref-anno ${gtf} --output alt_polya.bed --skip-feature-filter
    """
}

// splicing: frozen clusters (restricted file) quantified at every site + site-local novel clusters
process AGG_SPLICING {
    label 'process_low'

    input:
    path juncs, stageAs: 'junc/*'
    path restricted
    path gtf

    output:
    path 'splicing.bed', emit: bed
    path 'splicing.clusters_summary.tsv', emit: summary

    script:
    """
    ls \$PWD/junc/*.junc > juncfiles.txt
    xqtl_restricted.py extract ${restricted} restricted
    chroms=\$(cut -f1 ${gtf} | sort -u | paste -sd, -)
    xqtl_splicing_clusters.py frozen-refined --frozen restricted/splicing/frozen_clusters.tsv.gz --chroms "\$chroms" \\
        --out frozen_refined.txt --map frozen_map.tsv
    mkdir -p local frozen
    leafcutter_cluster_regtools_py3.py -j juncfiles.txt -r local -o leafcutter \\
        --maxintronlen 100000 --minclureads 30 --mincluratio 0.001
    leafcutter_cluster_regtools_py3.py -j juncfiles.txt -r frozen -o leafcutter -c frozen_refined.txt
    xqtl_splicing_clusters.py merge --frozen-numers frozen/leafcutter_perind_numers.counts.gz --frozen-map frozen_map.tsv \\
        --local-numers local/leafcutter_perind_numers.counts.gz --out splicing_numers.counts.gz \\
        --summary splicing.clusters_summary.tsv
    xqtl_assemble_bed.py splicing --input splicing_numers.counts.gz --ref-anno ${gtf} --output splicing.bed
    """
}

process AGG_INTRON_RETENTION {
    label 'process_low'

    input:
    path ir, stageAs: 'ir/*'
    path samples
    path gtf

    output:
    path 'intron_retention.bed', emit: bed

    script:
    """
    xqtl_assemble_bed.py intron-retention --samples ${samples} --input-dir ir --ref-anno ${gtf} --output intron_retention.bed
    """
}

process AGG_EDITING {
    label 'process_low'

    input:
    path levels, stageAs: 'levels/*'
    path samples
    tuple path(sites), path(site_genes)
    path gtf

    output:
    path 'rna_editing.bed', emit: bed
    path 'rna_editing.site_to_phenotype.tsv', emit: site_map

    script:
    // editing sites: coverage >= 20 reads in at least min(max(10, ceil(0.1 n)), n) samples; sites not clustered
    """
    n=\$(grep -c . ${samples})
    min_s=\$(awk -v n="\$n" 'BEGIN{m=int(0.1*n); if (m < 0.1*n) m++; if (m < 10) m=10; if (m > n) m=n; print m}')
    shared_samples_sites_matrix.py --path_to_edit_files levels --samples_file ${samples} \\
        --output_file edit_site_matrix.tsv --min_coverage 20 --min_samples "\$min_s"
    prepare_rna_editing_phenotypes.py --edit-matrix edit_site_matrix.tsv --site-to-gene ${site_genes} \\
        --output-matrix phenotype_matrix.tsv --output-site-map rna_editing.site_to_phenotype.tsv \\
        --correlation-threshold 0.9
    xqtl_assemble_bed.py rna-editing --input phenotype_matrix.tsv --ref-anno ${gtf} --output rna_editing.bed
    """
}

process AGG_STABILITY {
    label 'process_low'

    input:
    path counts, stageAs: 'fc/*'
    path samples
    path gtf

    output:
    path 'stability.bed', emit: bed

    script:
    """
    xqtl_assemble_bed.py stability --samples ${samples} --input-dir fc --ref-anno ${gtf} --output stability.bed --skip-feature-filter
    """
}

process QC_MERGE {
    label 'process_single'

    input:
    path rows, stageAs: 'rows/*'
    path quants, stageAs: 'salmon/*'
    path samples
    path gtf
    path gene_gc

    output:
    path 'qc_metrics.tsv', emit: qc

    script:
    """
    xqtl_picard_qc.py merge --rows rows/*.qcmetrics.tsv --salmon-dir salmon --samples ${samples} \\
        --gtf ${gtf} --gene-annot ${gene_gc} --output qc_metrics.tsv
    """
}

process MUSIC {
    label 'process_medium'

    input:
    path quants, stageAs: 'salmon/*'
    path samples
    path gtf
    path reference

    output:
    path 'cell_proportions.tsv', emit: props
    path 'cell_proportions_full.tsv', emit: full
    path 'maternal_flag.tsv', emit: flag

    script:
    """
    xqtl_bulk_tpm_matrix.py --gtf ${gtf} --salmon-dir salmon --samples ${samples} --output bulk_tpm_symbols.tsv
    xqtl_music.R --reference ${reference} --bulk bulk_tpm_symbols.tsv --cohort cohort --output-dir .
    mv cohort_cell_proportions_collapsed.tsv cell_proportions.tsv
    mv cohort_cell_proportions_full.tsv cell_proportions_full.tsv
    mv cohort_maternal_flag.tsv maternal_flag.tsv
    """
}

// sort, bgzip and index every BED; write the phenotypes_dir that --mode site maps
process PHENO_FINALIZE {
    label 'process_single'
    publishDir path: { "${params.outdir}" }, mode: params.publish_mode

    input:
    path tables, stageAs: 'in/*'

    output:
    path 'phenotypes', emit: dir

    script:
    """
    mkdir -p phenotypes
    for f in in/*; do
        b=\$(basename "\$f")
        case "\$b" in
          *.bed)
            { head -1 "\$f" | sed '1s/^chr\\t/#chr\\t/'; tail -n +2 "\$f" | sort -t \$'\\t' -k1,1 -k2,2n -k3,3n -k4,4; } \\
                | bgzip -c > "phenotypes/\$b.gz"
            tabix -p bed "phenotypes/\$b.gz" ;;
          *) cp -L "\$f" "phenotypes/\$b" ;;
        esac
    done
    {
      printf 'file\\tfeatures\\tsamples\\n'
      for f in phenotypes/*.bed.gz; do
        printf '%s\\t%s\\t%s\\n' "\$(basename "\$f")" "\$(gzip -dc "\$f" | tail -n +2 | wc -l)" "\$(gzip -dc "\$f" | head -1 | awk -F'\\t' '{print NF-4}')"
      done
    } > phenotypes/PHENOTYPES_SUMMARY.tsv
    """
}
