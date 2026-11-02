/*
 * --mode prepare_references : site indices built once from the public references embedded in the
 * image (refs_dir) and published to params.site_refs. Layout (read by the phenotyping workflow):
 *
 *   genome/genome.fa(.fai,.dict)          uncompressed genome for STAR, IRFinder-S, Picard, mpileup
 *   annotation/annotation.gtf             transcript annotation (header lines removed)
 *   annotation/exonic.gtf, intronic.gtf   constitutive exons / intronic gaps (featureCounts, stability)
 *   annotation/refFlat.txt                Picard CollectRnaSeqMetrics
 *   annotation/gene_gc_length.tsv         exonic GC content and length per gene (HCP subject bias)
 *   star/                                 STAR index (sjdbOverhang = rnaseq.read_length - 1)
 *   salmon/expression/                    Salmon index of the annotation transcripts (no decoys)
 *   salmon/txrevise.grp_{1,2}.{upstream,contained,downstream}/   Salmon indices of txrevise events
 *   irfinder/IRFinder/                    IRFinder-S reference (mapability exclusion included)
 *   editing/edit_sites.bed, edit_site_genes.tsv   (only when REDIportal sites are available)
 *   MANIFEST.tsv                          inputs (sha256), read length and tool versions
 */

process PREP_GENOME {
    label 'process_low'
    publishDir path: { "${params.site_refs}" }, mode: params.publish_mode

    input:
    path fasta

    output:
    tuple path('genome/genome.fa'), path('genome/genome.fa.fai'), path('genome/genome.dict'), emit: genome

    script:
    """
    mkdir -p genome
    if gzip -t ${fasta} 2>/dev/null; then gzip -dc ${fasta} > genome/genome.fa; else cat ${fasta} > genome/genome.fa; fi
    samtools faidx genome/genome.fa
    picard CreateSequenceDictionary R=genome/genome.fa O=genome/genome.dict
    """
}

process PREP_ANNOTATION {
    label 'process_low'
    publishDir path: { "${params.site_refs}" }, mode: params.publish_mode

    input:
    path gtf
    tuple path(fa), path(fai), path(dict)

    output:
    path 'annotation/annotation.gtf', emit: gtf
    tuple path('annotation/exonic.gtf'), path('annotation/intronic.gtf'), emit: exin
    path 'annotation/refFlat.txt', emit: refflat
    path 'annotation/gene_gc_length.tsv', emit: gene_gc

    script:
    // IRFinder-S's GTF parser reads '##gff-version' headers as a GFF switch, so headers are dropped
    """
    mkdir -p annotation
    if gzip -t ${gtf} 2>/dev/null; then gzip -dc ${gtf}; else cat ${gtf}; fi | grep -v '^#' > annotation/annotation.gtf
    exonic_intronic_from_gtf.py annotation/annotation.gtf annotation/exonic.gtf annotation/intronic.gtf --min-exon-fraction 0.8
    xqtl_picard_qc.py generate-refflat --gtf annotation/annotation.gtf --output annotation/refFlat.txt
    xqtl_picard_qc.py gene-gc-length --gtf annotation/annotation.gtf --fasta ${fa} --output annotation/gene_gc_length.tsv
    """
}

process STAR_INDEX {
    label 'process_high'
    publishDir path: { "${params.site_refs}" }, mode: params.publish_mode

    input:
    tuple path(fa), path(fai), path(dict)
    path gtf
    val read_length

    output:
    path 'star', emit: index

    script:
    """
    mkdir -p star
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir star \\
        --genomeFastaFiles ${fa} --sjdbGTFfile ${gtf} --sjdbOverhang ${(read_length as int) - 1} \\
        --genomeSAindexNbases ${params.star_sa_index_nbases}
    # genomeGenerate writes its Log.out into the current directory or the genome dir, depending on the STAR version
    if [ -f Log.out ]; then mv Log.out star/build.log; elif [ -f star/Log.out ]; then mv star/Log.out star/build.log; fi
    """
}

// one task per index: the annotation transcripts (expression) and the six txrevise event sets
process SALMON_INDEX {
    tag "${name}"
    label 'process_medium'
    publishDir path: { "${params.site_refs}/salmon" }, mode: params.publish_mode

    input:
    tuple val(name), path(models)
    tuple path(fa), path(fai), path(dict)

    output:
    tuple val(name), path("${name}"), emit: index

    script:
    """
    if gzip -t ${models} 2>/dev/null; then gzip -dc ${models} > models.txt; else cp ${models} models.txt; fi
    gffread models.txt -g ${fa} -w transcripts.fa
    test -s transcripts.fa
    salmon index -t transcripts.fa -i ${name} -p ${task.cpus}
    rm -f transcripts.fa models.txt
    """
}

process IRFINDER_REF {
    label 'process_high'
    publishDir path: { "${params.site_refs}/irfinder" }, mode: params.publish_mode

    input:
    tuple path(fa), path(fai), path(dict)
    path gtf
    path star
    val read_length

    output:
    path 'IRFinder', emit: ref

    script:
    // only REF/IRFinder is needed by `IRFinder BAM`; STAR is linked (-l), not copied
    """
    IRFinder BuildRefFromSTARRef -r REF -x ${star} -f ${fa} -g ${gtf} -l -t ${task.cpus} -n ${read_length}
    mv REF/IRFinder IRFinder
    cp REF/*.log IRFinder/ 2>/dev/null || true
    """
}

process EDIT_SITES {
    label 'process_single'
    publishDir path: { "${params.site_refs}" }, mode: params.publish_mode

    input:
    path sites
    path gtf

    output:
    tuple path('editing/edit_sites.bed'), path('editing/edit_site_genes.tsv'), emit: sites

    script:
    // sites are restricted to chromosomes present in the annotation (keeps chr22 tests small)
    """
    mkdir -p editing
    cut -f1 ${gtf} | sort -u > chroms.txt
    if gzip -t ${sites} 2>/dev/null; then gzip -dc ${sites}; else cat ${sites}; fi \\
        | awk 'BEGIN{FS=OFS="\\t"} NR==FNR{c[\$1]=1; next} (\$1 in c)' chroms.txt - > editing/edit_sites.bed
    test -s editing/edit_sites.bed
    map_edit_sites_to_genes.py editing/edit_sites.bed ${gtf} editing/edit_site_genes.tsv
    """
}

process REFS_MANIFEST {
    label 'process_single'
    publishDir path: { "${params.site_refs}" }, mode: params.publish_mode

    input:
    path inputs, stageAs: 'inputs/*'
    val read_length
    val done

    output:
    path 'MANIFEST.tsv', emit: manifest

    script:
    """
    {
      printf 'key\\tvalue\\n'
      printf 'pipeline_version\\t%s\\n' '${workflow.manifest.version}'
      printf 'read_length\\t%s\\n' '${read_length}'
      printf 'star_sjdb_overhang\\t%s\\n' '${(read_length as int) - 1}'
      printf 'built\\t%s\\n' "\$(date -u +%Y-%m-%dT%H:%M:%SZ)"
      for f in inputs/*; do printf 'sha256:%s\\t%s\\n' "\$(basename "\$f")" "\$(sha256sum "\$f" | cut -d' ' -f1)"; done
      xqtl_tool_versions.sh 2>/dev/null | awk 'BEGIN{FS=OFS="\\t"} NR>1 && NF>=2 {print "version:"\$1, \$2}'
    } > MANIFEST.tsv
    """
}
