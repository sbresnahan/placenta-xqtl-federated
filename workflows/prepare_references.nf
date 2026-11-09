/*
 * --mode prepare_references : build the site indices (STAR, Salmon, IRFinder-S, Picard and
 * stability annotations) from the public references and publish them to params.site_refs.
 * Run once per site and read length; the phenotyping workflow checks MANIFEST.tsv.
 */

include { PREP_GENOME ; PREP_ANNOTATION ; STAR_INDEX ; SALMON_INDEX ; IRFINDER_REF ; EDIT_SITES ;
          REFS_MANIFEST } from '../modules/local/references.nf'

// reference inputs: explicit params win, else the files embedded in the image under refs_dir
def refInput(c, key, rel) {
    return c[key] ?: "${c.refs_dir}/${rel}".toString()
}

// Only the upstream / downstream event sets feed phenotypes (alternative TSS / polyA), as in
// Pantry; the 'contained' sets are distributed with the references but not indexed.
def txreviseSets() {
    return ['grp_1', 'grp_2'].collectMany { g -> ['upstream', 'downstream'].collect { p -> "txrevise.${g}.${p}".toString() } }
}

workflow PREPARE_REFERENCES {
    take:
    cfg     // value channel: resolved cohort parameters (map)

    main:
    fasta = cfg.map { c -> file(refInput(c, 'genome_fasta', 'genome/GRCh38_no_alt_analysis_set.fa.gz'), checkIfExists: true) }
    gtf_in = cfg.map { c -> file(refInput(c, 'annotation_gtf', 'annotation/HPLRv2.0.annotated.PANTRY.gtf.gz'), checkIfExists: true) }
    read_length = cfg.map { c -> c.rnaseq.read_length as int }

    PREP_GENOME(fasta)
    genome = PREP_GENOME.out.genome
    PREP_ANNOTATION(gtf_in, genome)
    gtf = PREP_ANNOTATION.out.gtf

    STAR_INDEX(genome, gtf, read_length)

    tx_models = cfg.flatMap { c ->
        def dir = refInput(c, 'txrevise_dir', 'txrevise')
        txreviseSets().collect { s -> tuple(s, file("${dir}/${s}.gff3.gz", checkIfExists: true)) }
    }
    SALMON_INDEX(gtf.map { g -> tuple('expression', g) }.mix(tx_models), genome)

    IRFINDER_REF(genome, gtf, STAR_INDEX.out.index, read_length)

    sites = cfg.flatMap { c ->
        def f = file(refInput(c, 'editing_sites', 'editing/rediportal_hg38.bed.gz'))
        if (!f.exists()) {
            log.warn("No REDIportal editing sites at ${f}: the rna_editing modality will be skipped (docs/REFERENCES.md)")
        }
        f.exists() ? [f] : []
    }
    EDIT_SITES(sites, gtf)

    done = STAR_INDEX.out.index.mix(SALMON_INDEX.out.index, IRFINDER_REF.out.ref, PREP_ANNOTATION.out.refflat,
        PREP_ANNOTATION.out.gene_gc, EDIT_SITES.out.sites).collect().map { _l -> true }
    inputs = fasta.mix(gtf_in, tx_models.map { t -> t[1] }, sites).collect()
    REFS_MANIFEST(inputs, read_length, done)

    emit:
    manifest = REFS_MANIFEST.out.manifest
}
