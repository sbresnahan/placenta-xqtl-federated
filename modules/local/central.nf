/*
 * --mode central (MDACC): meta-analysis of site exports. Scripts: bin/xqtl_central.py (ingest, meta,
 * features, fine-mapping/colocalization inputs, G x C, G x E, TWAS inputs, summary), bin/xqtl_finemap.py
 * (SuShiE), bin/xqtl_coloc.R (SuSiE-coloc, ColocBoost), bin/xqtl_fusion_wgt.R (FUSION weights).
 * Study export folders are staged in studies.tsv order as study01/<dir>, study02/<dir>, ...; the staged
 * relative paths recorded by CEN_INGEST in contributions.tsv are valid in every later task.
 * Outputs are published to ${outdir}/central/.
 */

process CEN_INGEST {
    label 'process_medium'
    publishDir path: "${params.outdir}/central", mode: params.publish_mode

    input:
    path studies
    path exports, stageAs: 'study??/*'
    path config

    output:
    path 'ingest', emit: ingest
    path 'ingest/units.tsv', emit: units

    script:
    """
    xqtl_central.py ingest --studies ${studies} --config ${config} --staged-prefix study --out-dir ingest
    """
}

process CEN_META {
    tag "${modality} ${chrom}"
    label 'process_medium'

    input:
    tuple val(modality), val(chrom)
    path ingest
    path exports, stageAs: 'study??/*'

    output:
    path "meta/${modality}_${chrom}.parquet", emit: meta

    script:
    """
    mkdir -p meta
    xqtl_central.py meta-nominal --contributions ${ingest}/contributions.tsv --harmonization ${ingest}/harmonization \\
        --modality ${modality} --chrom ${chrom} --out meta/${modality}_${chrom}.parquet
    """
}

process CEN_FEATURES {
    label 'process_medium'
    publishDir path: "${params.outdir}/central", mode: params.publish_mode, pattern: '{features.tsv.gz,meta}'

    input:
    path meta, stageAs: 'meta/*'
    path ingest
    path exports, stageAs: 'study??/*'
    path config

    output:
    path 'features.tsv.gz', emit: features
    path 'meta', emit: meta

    script:
    """
    xqtl_central.py features --contributions ${ingest}/contributions.tsv --harmonization ${ingest}/harmonization \\
        --meta-dir meta --config ${config} --out features.tsv.gz
    """
}

process CEN_FINEMAP_PREP {
    label 'process_medium'
    publishDir path: "${params.outdir}/central/finemap", mode: params.publish_mode, pattern: 'loci.tsv'

    input:
    path meta, stageAs: 'meta/*'
    path features
    path ingest
    path exports, stageAs: 'study??/*'
    path config

    output:
    path 'loci.tsv', emit: loci
    path 'loci/*.npz', optional: true, emit: npz

    script:
    """
    xqtl_central.py finemap-prep --contributions ${ingest}/contributions.tsv --harmonization ${ingest}/harmonization \\
        --meta-dir meta --features ${features} --config ${config} --out-dir loci
    mv loci/loci.tsv loci.tsv
    """
}

process CEN_FINEMAP {
    tag "shard ${shard}"
    label 'process_medium'

    input:
    tuple val(shard), path(npz, stageAs: 'loci/*')
    val f

    output:
    path "fm_${shard}.*.tsv", emit: tables

    script:
    """
    xqtl_finemap.py --loci loci/*.npz --L ${f.L} --purity ${f.purity} --out-prefix fm_${shard}
    """
}

process CEN_CONCAT {
    tag "${stage}"
    label 'process_single'
    publishDir path: { "${params.outdir}/central/${stage}" }, mode: params.publish_mode

    input:
    tuple val(stage), val(kind), path(parts, stageAs: 'parts/*')

    output:
    tuple val(kind), path("${kind}.tsv"), emit: table

    script:
    """
    xqtl_central.py concat --out ${kind}.tsv parts/*
    """
}

process CEN_COLOC_PREP {
    label 'process_medium'
    publishDir path: "${params.outdir}/central/coloc", mode: params.publish_mode, pattern: 'coloc_tasks.tsv'

    input:
    path catalog
    path gwas, stageAs: 'gwas/*'
    path loci
    path npz, stageAs: 'loci/*'
    path config

    output:
    path 'coloc_tasks.tsv', emit: tasks
    path 'coloc_in', emit: dir

    script:
    """
    xqtl_central.py coloc-prep --catalog ${catalog} --loci ${loci} --loci-dir loci --config ${config} \\
        --out-dir coloc_in --staged-dir gwas
    cp coloc_in/coloc_tasks.tsv coloc_tasks.tsv
    """
}

process CEN_COLOC {
    tag "shard ${shard}"
    label 'process_medium'

    input:
    tuple val(shard), val(task_ids)
    path dir
    val c

    output:
    path "co_${shard}.coloc.tsv", emit: coloc
    path "co_${shard}.coloc_pairs.tsv", emit: pairs
    path "co_${shard}.colocboost.tsv", optional: true, emit: colocboost

    script:
    def cb = c.colocboost ? "--colocboost --M ${c.M}" : ''
    """
    printf '%s\\n' ${task_ids.join(' ')} | paste -sd, - > ids.txt
    xqtl_coloc.R --tasks ${dir}/coloc_tasks.tsv --task-dir ${dir} --task-ids "\$(cat ids.txt)" \\
        --pp-h4 ${c.pp_h4} --L ${c.L} ${cb} --out-prefix co_${shard}
    """
}

process CEN_INTERACTIONS {
    label 'process_medium'
    publishDir path: "${params.outdir}/central", mode: params.publish_mode

    input:
    path ingest
    path exports, stageAs: 'study??/*'
    path config

    output:
    path 'interactions', emit: dir

    script:
    """
    for s in gxc gxe-tier1 gxe-tier2; do
        xqtl_central.py \$s --contributions ${ingest}/contributions.tsv --harmonization ${ingest}/harmonization \\
            --config ${config} --out-dir interactions
    done
    """
}

process CEN_TWAS {
    label 'process_low'
    publishDir path: "${params.outdir}/central", mode: params.publish_mode

    input:
    path meta, stageAs: 'meta/*'
    path features
    path weights
    path loci
    path ingest
    path exports, stageAs: 'study??/*'
    path config

    output:
    path 'twas', emit: dir

    script:
    def w = weights && loci ? "--weights ${weights} --loci ${loci}" : ''
    """
    xqtl_central.py twas-inputs --contributions ${ingest}/contributions.tsv --meta-dir meta --features ${features} \\
        --config ${config} --out-dir twas ${w}
    xqtl_fusion_wgt.R --twas-dir twas
    """
}

process CEN_SUMMARY {
    label 'process_low'
    publishDir path: "${params.outdir}/central", mode: params.publish_mode

    input:
    path features
    path loci
    path fm, stageAs: 'fm/*'
    path co, stageAs: 'co/*'
    path interactions
    path annotations
    path beds, stageAs: 'annot/*'
    path ingest
    path exports, stageAs: 'study??/*'
    path config

    output:
    path 'summary', emit: dir

    script:
    def ann = annotations ? "--annotations ${annotations} --staged-dir annot" : ''
    def lo = loci ? "--loci ${loci}" : ''
    """
    args=""
    for k in pips cs status; do if [ -s fm/\$k.tsv ]; then args="\$args --\$k fm/\$k.tsv"; fi; done
    for k in coloc colocboost; do if [ -s co/\$k.tsv ]; then args="\$args --\$k co/\$k.tsv"; fi; done
    # shellcheck disable=SC2086
    xqtl_central.py summary --contributions ${ingest}/contributions.tsv --config ${config} --out-dir summary \\
        --features ${features} ${lo} \$args ${ann} \\
        --gxc ${interactions}/gxc_annotation.tsv.gz --tier1 ${interactions}/gxe_tier1_meta.tsv.gz
    cp ${ingest}/ingest_report.tsv ${ingest}/harmonization_summary.tsv summary/
    """
}
