/*
 * --mode export : rebuild the shareable export (A-F, I, J) from a finished site run, e.g. after a
 * change of export.format / export.round_digits or when the leak scan needs to be re-run.
 * Reads ${outdir}/site/: ancestry/cohort_strata.tsv, <ANC>/ staging files, genotype_qc/
 * cohort.vcf_samples.txt and restricted_version.txt. Nothing is recomputed.
 */

include { EXPORT_STRATUM ; EXPORT_FINALIZE } from '../modules/local/export.nf'

def mappedStrata(strataFile) {
    def lines = strataFile.readLines()
    def hdr = lines[0].tokenize('\t')
    return lines.drop(1).collect { l -> l.tokenize('\t') }
        .findAll { r -> r[hdr.indexOf('status')] == 'mapped' }
        .collect { r -> r[hdr.indexOf('ancestry')] }
}

workflow EXPORT_ONLY {
    take:
    cfg     // value channel: resolved cohort parameters (map)

    main:
    site = cfg.map { c -> file("${c.outdir}/site", checkIfExists: true) }
    strata = site.map { s -> file("${s}/ancestry/cohort_strata.tsv", checkIfExists: true) }
    staging = strata.combine(site).flatMap { st, s ->
        mappedStrata(st).collect { anc ->
            def staged = files("${s}/${anc}/${anc}*", type: 'any')
            if (!staged) {
                error("--mode export: no site outputs for stratum ${anc} in ${s}/${anc}")
            }
            tuple(anc, staged)
        }
    }
    phen = cfg.map { c -> c.phenotypes_dir ?: "${c.outdir}/phenotypes" }
    qc_metrics = phen.map { d -> file("${d}/qc_metrics.tsv", checkIfExists: true) }
    cell_props = phen.map { d -> file("${d}/cell_proportions.tsv").exists() ? file("${d}/cell_proportions.tsv") : [] }
    ex = cfg.combine(site).map { c, s ->
        c.export + [cohort_id: c.cohort_id, image_sha256: c.image_sha256 ?: 'unknown',
                    restricted_version: file("${s}/restricted_version.txt", checkIfExists: true).text.trim(),
                    pipeline_version: workflow.manifest.version]
    }
        .first()
    EXPORT_STRATUM(staging, qc_metrics, cell_props, ex)
    ids = cfg.combine(site).map { c, s ->
        [file(c.metadata, checkIfExists: true), file(c.samplesheet, checkIfExists: true),
         file("${s}/genotype_qc/cohort.vcf_samples.txt", checkIfExists: true)]
    }
        .first()
    EXPORT_FINALIZE(EXPORT_STRATUM.out.dir.collect(), strata, ids, ex)

    emit:
    export = EXPORT_FINALIZE.out.export
}
