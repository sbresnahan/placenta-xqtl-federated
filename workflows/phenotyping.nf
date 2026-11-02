/*
 * Phenotyping: FASTQs -> the nine phenotype BEDs, count BEDs, qc_metrics.tsv and
 * cell_proportions.tsv (the phenotypes_dir contract in docs/EXPORT_SPEC.md), published to
 * ${outdir}/phenotypes/. Used by --mode site (followed by mapping) and --mode phenotyping_only
 * (followed by normalization over all samples, no export).
 *
 * Indices come from params.site_refs (--mode prepare_references). Splicing is quantified on the
 * frozen leafCutter clusters from the restricted file plus site-local novel clusters.
 */

include { STAR_ALIGN ; SALMON_EXPRESSION ; SALMON_TXREVISE ; REGTOOLS_JUNCTIONS ; FEATURECOUNTS ;
          EDITING_LEVELS ; IRFINDER ; PICARD_QC ; AGG_EXPRESSION ; AGG_ALT_ENDS ; AGG_SPLICING ;
          AGG_INTRON_RETENTION ; AGG_EDITING ; AGG_STABILITY ; QC_MERGE ; MUSIC ;
          PHENO_FINALIZE } from '../modules/local/phenotyping.nf'
include { NORMALIZE as NORMALIZE_ALL ; NORMALIZE as NORMALIZE_ALL_ISOFORM } from '../modules/local/mapping.nf'

def siteRef(c, rel) {
    return file("${c.site_refs}/${rel}", checkIfExists: true)
}

def hasEditing(c) {
    return file("${c.site_refs}/editing/edit_sites.bed").exists()
}

def musicReference(c) {
    return file(c.music_reference ?: "${c.refs_dir}/music/placenta_music_reference.rds")
}

// modalities produced at this site: cohort.yml mapping.modalities, minus those a 3' library
// cannot measure and minus rna_editing when no REDIportal sites were prepared
def phenoModalities(c) {
    def m = c.mapping.modalities as List
    if (c.rnaseq.library_type == '3prime') {
        m = m.findAll { x -> x in ['expression', 'stability'] }
    }
    if (!hasEditing(c)) {
        m = m.findAll { x -> x != 'rna_editing' }
    }
    return (['expression'] + m).unique().join(',')
}

// keep the items of `ch` only when modality `name` is produced (mods: value channel, comma string)
def onlyIf(ch, mods, name) {
    return ch.combine(mods).filter { t -> name in t[-1].tokenize(',') }.map { t -> t.size() == 2 ? t[0] : t[0..-2] }
}

// FASTQ paths in the samplesheet may be relative to the samplesheet's directory
def fastqPath(base, p) {
    return p.startsWith('/') ? file(p, checkIfExists: true) : file(base.resolve(p).toString(), checkIfExists: true)
}

def readsOf(c) {
    def base = file(c.samplesheet).parent
    return file(c.samplesheet).splitCsv(header: true).collect { r ->
        tuple(r.sample_id, r.fastq_2 ? [fastqPath(base, r.fastq_1), fastqPath(base, r.fastq_2)] : [fastqPath(base, r.fastq_1)])
    }
}

def isIsoform(m) {
    return m in ['isoform_expression', 'isoform_usage']
}

workflow PHENOTYPING {
    take:
    cfg     // value channel: resolved cohort parameters (map)

    main:
    mods = cfg.map { c -> phenoModalities(c) }
    mods.view { m -> "Phenotyping modalities: ${m}" }
    rnaseq = cfg.map { c -> c.rnaseq }
    restricted = cfg.map { c -> file(c.restricted, checkIfExists: true) }

    genome = cfg.map { c -> tuple(siteRef(c, 'genome/genome.fa'), siteRef(c, 'genome/genome.fa.fai'), siteRef(c, 'genome/genome.dict')) }
    gtf = cfg.map { c -> siteRef(c, 'annotation/annotation.gtf') }
    exin = cfg.map { c -> tuple(siteRef(c, 'annotation/exonic.gtf'), siteRef(c, 'annotation/intronic.gtf')) }
    refflat = cfg.map { c -> siteRef(c, 'annotation/refFlat.txt') }
    gene_gc = cfg.map { c -> siteRef(c, 'annotation/gene_gc_length.tsv') }
    star = cfg.map { c -> siteRef(c, 'star') }
    salmon_expr = cfg.map { c -> siteRef(c, 'salmon/expression') }
    salmon_tx = cfg.map { c ->
        ['grp_1.upstream', 'grp_2.upstream', 'grp_1.downstream', 'grp_2.downstream'].collect { s -> siteRef(c, "salmon/txrevise.${s}") }
    }
    irf_ref = cfg.map { c -> siteRef(c, 'irfinder/IRFinder') }
    edit_sites = cfg.flatMap { c -> hasEditing(c) ? [tuple(siteRef(c, 'editing/edit_sites.bed'), siteRef(c, 'editing/edit_site_genes.tsv'))] : [] }

    reads = cfg.flatMap { c -> readsOf(c) }
    samples = reads.map { s, _r -> s }.collectFile(name: 'samples.txt', newLine: true, sort: true)

    // ---------------- per sample ----------------
    STAR_ALIGN(reads, star)
    bam = STAR_ALIGN.out.bam
    SALMON_EXPRESSION(reads, salmon_expr)
    alt_reads = reads.combine(mods).filter { t -> t[-1].tokenize(',').any { m -> m in ['alt_tss', 'alt_polya'] } }
        .map { s, r, _m -> tuple(s, r) }
    SALMON_TXREVISE(alt_reads, salmon_tx)
    REGTOOLS_JUNCTIONS(onlyIf(bam, mods, 'splicing'))
    FEATURECOUNTS(onlyIf(bam, mods, 'stability'), exin, rnaseq)
    EDITING_LEVELS(onlyIf(bam, mods, 'rna_editing'), genome, edit_sites.first())
    IRFINDER(onlyIf(bam, mods, 'intron_retention'), irf_ref)
    PICARD_QC(bam, genome, refflat, rnaseq)

    // ---------------- cohort level ----------------
    quants = SALMON_EXPRESSION.out.quant.map { _s, d -> d }.collect()
    AGG_EXPRESSION(quants, samples, gtf)
    AGG_ALT_ENDS(SALMON_TXREVISE.out.quant.map { _s, d -> d }.collect(), samples, gtf)
    AGG_SPLICING(REGTOOLS_JUNCTIONS.out.junc.collect(), restricted, gtf)
    AGG_INTRON_RETENTION(IRFINDER.out.ir.collect(), samples, gtf)
    AGG_EDITING(EDITING_LEVELS.out.levels.collect(), samples, edit_sites.first(), gtf)
    AGG_STABILITY(FEATURECOUNTS.out.counts.flatten().collect(), samples, gtf)
    QC_MERGE(PICARD_QC.out.row.collect(), quants, samples, gtf, gene_gc)
    music_ref = cfg.flatMap { c ->
        def f = musicReference(c)
        if (!f.exists()) {
            log.warn("MuSiC reference not found (${f}): cell_proportions.tsv will not be written")
        }
        f.exists() ? [f] : []
    }
    MUSIC(quants, samples, gtf, music_ref)

    tables = AGG_EXPRESSION.out.expression
        .mix(AGG_EXPRESSION.out.counts.flatten().filter { f -> f.name.startsWith('expression.') })
        .mix(onlyIf(AGG_EXPRESSION.out.isoform_expression, mods, 'isoform_expression'))
        .mix(onlyIf(AGG_EXPRESSION.out.counts.flatten().filter { f -> f.name.startsWith('isoform_expression.') }, mods, 'isoform_expression'))
        .mix(onlyIf(AGG_EXPRESSION.out.isoform_usage, mods, 'isoform_usage'))
        .mix(onlyIf(AGG_ALT_ENDS.out.alt_tss, mods, 'alt_tss'))
        .mix(onlyIf(AGG_ALT_ENDS.out.alt_polya, mods, 'alt_polya'))
        .mix(AGG_SPLICING.out.bed, AGG_SPLICING.out.summary)
        .mix(AGG_INTRON_RETENTION.out.bed)
        .mix(AGG_EDITING.out.bed, AGG_EDITING.out.site_map)
        .mix(AGG_STABILITY.out.bed)
        .mix(QC_MERGE.out.qc)
        .mix(MUSIC.out.props, MUSIC.out.full, MUSIC.out.flag)
        .collect()
    PHENO_FINALIZE(tables)

    emit:
    dir = PHENO_FINALIZE.out.dir
}

// --mode phenotyping_only: phenotype BEDs plus normalized matrices over all samples (no strata,
// no genotypes, no export). Normalized BEDs are published to ${outdir}/normalized/.
workflow PHENOTYPING_ONLY {
    take:
    cfg

    main:
    PHENOTYPING(cfg)
    metadata = cfg.map { c -> file(c.metadata, checkIfExists: true) }
    norm_in = PHENOTYPING.out.dir.combine(cfg).flatMap { d, c ->
        c.mapping.modalities.findAll { m -> d.resolve("${m}.bed.gz").exists() }.collect { m ->
            def counts = d.resolve("${m}.counts.bed.gz")
            tuple('ALL', m, d.resolve("${m}.bed.gz"), counts.exists() ? counts : [])
        }
    }
        .combine(metadata)
    NORMALIZE_ALL(norm_in.filter { t -> !isIsoform(t[1]) }.map { anc, m, bed, counts, smp -> tuple(anc, m, bed, counts, [], smp) })
    expr_outliers = NORMALIZE_ALL.out.outliers.filter { t -> t[1] == 'expression' }.map { anc, _m, o -> tuple(anc, o) }
    NORMALIZE_ALL_ISOFORM(norm_in.filter { t -> isIsoform(t[1]) }.combine(expr_outliers, by: 0)
        .map { anc, m, bed, counts, smp, o -> tuple(anc, m, bed, counts, o, smp) })

    emit:
    dir = PHENOTYPING.out.dir
    normalized = NORMALIZE_ALL.out.bed.mix(NORMALIZE_ALL_ISOFORM.out.bed).map { _a, _m, bed, _g -> bed }.collect()
}
