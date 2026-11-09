/*
 * --mode imputation_prep : genotypes.raw (VCF/BCF, PLINK 1 bed/bim/fam or PLINK 2 pgen/pvar/psam;
 * array, WES or WGS) -> ${outdir}/imputation_prep/chr1..22.vcf.gz + report + TOPMED_CHECKLIST.md.
 * GRCh37 input is lifted to GRCh38 first. The upload to the TOPMed Imputation Server is manual.
 */

include { IP_GENOME ; IP_CONVERT ; IP_LIFTOVER ; IP_FILTER ; IP_REPORT ; IP_PUBLISH } from '../modules/local/imputation_prep.nf'

// [kind, files, prefix (staged name passed to the script), build]
def rawInput(c) {
    def r = c.genotypes.raw.toString()
    def build = c.genotypes.raw_build ?: 'GRCh38'
    if (r ==~ /.*\.(vcf|vcf\.gz|vcf\.bgz|bcf)$/) {
        def f = file(r, checkIfExists: true)
        return ['vcf', [f], f.name, build]
    }
    if (file("${r}.bed").exists()) {
        return ['bed', ['bed', 'bim', 'fam'].collect { e -> file("${r}.${e}", checkIfExists: true) }, file(r).name, build]
    }
    if (file("${r}.pgen").exists()) {
        def pvar = file("${r}.pvar").exists() ? file("${r}.pvar") : file("${r}.pvar.zst", checkIfExists: true)
        return ['pgen', [file("${r}.pgen"), pvar, file("${r}.psam", checkIfExists: true)], file(r).name, build]
    }
    error("genotypes.raw: ${r} is neither a VCF/BCF file nor a PLINK prefix (.bed/.bim/.fam or .pgen/.pvar/.psam)")
}

def autosomes(chromsFile) {
    return chromsFile.readLines().drop(1).collect { l -> l.tokenize('\t')[0] }
        .findAll { ch -> ch ==~ /chr([1-9]|1[0-9]|2[0-2])/ }
}

workflow IMPUTATION_PREP {
    take:
    cfg     // value channel: resolved cohort parameters (map)

    main:
    fasta = cfg.map { c -> file(c.genome_fasta ?: "${c.refs_dir}/genome/GRCh38_no_alt_analysis_set.fa.gz", checkIfExists: true) }
    fasta_idx = fasta.map { f -> ['fai', 'gzi'].collect { e -> file("${f}.${e}") }.findAll { x -> x.exists() } }
    build = cfg.map { c -> c.genotypes.raw_build ?: 'GRCh38' }
    IP_GENOME(fasta, fasta_idx, build)

    IP_CONVERT(cfg.map { c -> rawInput(c) })

    lifted = cfg.filter { c -> c.genotypes.raw_build == 'GRCh37' }
    chain = lifted.map { c -> file(c.liftover_chain ?: "${c.refs_dir}/liftover/hg19ToHg38.over.chain.gz", checkIfExists: true) }
    IP_LIFTOVER(IP_CONVERT.out.vcf.combine(lifted).map { v, i, _c -> tuple(v, i) }, IP_GENOME.out.genome, chain)

    // GRCh38 input goes straight to filtering; GRCh37 input after liftover
    hg38 = IP_CONVERT.out.vcf.combine(IP_CONVERT.out.chroms).combine(cfg)
        .filter { _v, _i, _ch, c -> c.genotypes.raw_build != 'GRCh37' }
        .map { v, i, ch, _c -> tuple(v, i, ch) }
    to_filter = hg38.mix(IP_LIFTOVER.out.vcf.combine(IP_LIFTOVER.out.chroms))
        .flatMap { v, i, ch -> autosomes(ch).collect { chrom -> tuple(chrom, v, i) } }
    ip = cfg.map { c -> c.imputation_prep }
    IP_FILTER(to_filter, IP_GENOME.out.genome, ip)

    metadata = cfg.map { c -> c.metadata ? file(c.metadata, checkIfExists: true) : [] }
    lift_counts = IP_LIFTOVER.out.counts.ifEmpty([])
    IP_REPORT(IP_FILTER.out.counts.collect(), IP_CONVERT.out.smiss, IP_CONVERT.out.sexcheck, lift_counts, metadata, cfg)
    IP_PUBLISH(IP_FILTER.out.vcf.flatMap { _ch, v, t -> [v, t] }.collect(), IP_REPORT.out.report)

    emit:
    vcfs = IP_PUBLISH.out.vcfs
    checklist = IP_REPORT.out.checklist
}
