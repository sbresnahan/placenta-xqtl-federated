/*
 * --mode central (MDACC): site exports listed in central.studies -> ${outdir}/central/
 *   ingest (checksums, versions, schemas, allele harmonization) -> nominal IVW meta per modality x chromosome
 *   -> feature-level ACAT + q-values -> SuShiE multi-ancestry fine-mapping -> colocalization with the GWAS
 *   catalog (optional) ; G x C and G x E synthesis ; OTTERS / FUSION TWAS inputs ; summary + PIP enrichment.
 * See docs/CENTRAL_ANALYSIS.md.
 */

include { CEN_INGEST ; CEN_META ; CEN_FEATURES ; CEN_FINEMAP_PREP ; CEN_FINEMAP ; CEN_CONCAT as CEN_CONCAT_FM ;
          CEN_CONCAT as CEN_CONCAT_CO ; CEN_COLOC_PREP ; CEN_COLOC ; CEN_INTERACTIONS ; CEN_TWAS ;
          CEN_SUMMARY } from '../modules/local/central.nf'

// rows of a TSV with a header, relative paths resolved against the table's folder
def tsvRows(path) {
    def f = file(path, checkIfExists: true)
    def lines = f.readLines().findAll { l -> l.trim() }
    def head = lines[0].tokenize('\t')*.trim()
    return lines.drop(1).collect { l ->
        def v = l.split('\t', -1)*.trim()
        def m = [:]
        head.eachWithIndex { h, i -> m[h] = i < v.size() ? v[i] : '' }
        m['_base'] = f.parent
        m
    }
}

def resolvePath(String p, base) {
    return p.startsWith('/') ? file(p, checkIfExists: true) : file("${base}/${p}", checkIfExists: true)
}

workflow CENTRAL {
    take:
    cfg       // value channel: resolved parameters (map)
    config    // value channel: resolved_params.json (read by the scripts)

    main:
    studies = cfg.map { c -> file(c.central.studies, checkIfExists: true) }
    exports = cfg.map { c -> tsvRows(c.central.studies).collect { r -> resolvePath(r.export_dir, r._base) } }

    CEN_INGEST(studies, exports, config)
    ingest = CEN_INGEST.out.ingest

    units = CEN_INGEST.out.units.flatMap { f ->
        f.readLines().drop(1).findAll { l -> l.trim() }.collect { l -> tuple(l.tokenize('\t')[0], l.tokenize('\t')[1]) }
    }
    CEN_META(units, ingest, exports)
    meta = CEN_META.out.meta.collect()

    CEN_FEATURES(meta, ingest, exports, config)
    features = CEN_FEATURES.out.features

    // ---- fine-mapping -------------------------------------------------------------------------------
    fm_on = cfg.filter { c -> c.central.finemap.enabled }
    // a list-valued item would be flattened by combine: wrap it
    fm_meta = fm_on.map { _c -> 1 }.combine(meta.map { m -> [m] }).map { _x, m -> m }
    CEN_FINEMAP_PREP(fm_meta, features, ingest, exports, config)
    shard_size = cfg.map { c -> c.central.finemap.shard_size as int }
    npz_sorted = CEN_FINEMAP_PREP.out.npz.flatten().toSortedList { a, b -> a.name <=> b.name }
    shards = shard_size.combine(npz_sorted.map { l -> [l] })
        .flatMap { n, files -> files.collate(n).withIndex().collect { chunk, i -> tuple(i + 1, chunk) } }
    CEN_FINEMAP(shards, cfg.map { c -> c.central.finemap })
    fm_parts = CEN_FINEMAP.out.tables.flatten()
        .map { f -> tuple(f.name.tokenize('.')[-2], f) }
        .groupTuple()
        .map { kind, files -> tuple('finemap', kind, files) }
    CEN_CONCAT_FM(fm_parts)
    fm_tables = CEN_CONCAT_FM.out.table.map { _k, f -> f }.collect().ifEmpty([])
    weights = CEN_CONCAT_FM.out.table.filter { k, _f -> k == 'weights' }.map { _k, f -> f }.ifEmpty([])
    loci = CEN_FINEMAP_PREP.out.loci.ifEmpty([])

    // ---- colocalization (needs central.gwas_catalog) ---------------------------------------------
    co_cfg = cfg.filter { c -> c.central.gwas_catalog && c.central.coloc.enabled && c.central.finemap.enabled }
    catalog = co_cfg.map { c -> file(c.central.gwas_catalog, checkIfExists: true) }
    gwas = co_cfg.map { c -> tsvRows(c.central.gwas_catalog).collect { r -> resolvePath(r.path, r._base) } }
    CEN_COLOC_PREP(catalog, gwas, CEN_FINEMAP_PREP.out.loci, CEN_FINEMAP_PREP.out.npz.collect().ifEmpty([]), config)
    co_shards = CEN_COLOC_PREP.out.tasks.combine(co_cfg)
        .flatMap { f, c ->
            def ids = f.readLines().drop(1).collect { l -> l.tokenize('\t') }.findAll { x -> x[-1] == 'prepared' }.collect { x -> x[0] }
            ids.collate(c.central.coloc.shard_size as int).withIndex().collect { chunk, i -> tuple(i + 1, chunk) }
        }
    CEN_COLOC(co_shards, CEN_COLOC_PREP.out.dir.first(), co_cfg.map { c -> c.central.coloc })
    co_parts = CEN_COLOC.out.coloc.map { f -> tuple('coloc', f) }
        .mix(CEN_COLOC.out.pairs.map { f -> tuple('coloc_pairs', f) })
        .mix(CEN_COLOC.out.colocboost.map { f -> tuple('colocboost', f) })
        .groupTuple()
        .map { kind, files -> tuple('coloc', kind, files) }
    CEN_CONCAT_CO(co_parts)
    co_tables = CEN_CONCAT_CO.out.table.map { _k, f -> f }.collect().ifEmpty([])

    // ---- G x C, G x E, TWAS inputs, summary ------------------------------------------------------
    CEN_INTERACTIONS(ingest, exports, config)
    CEN_TWAS(meta, features, weights, loci, ingest, exports, config)

    annotations = cfg.map { c -> c.central.annotations ? file(c.central.annotations, checkIfExists: true) : [] }
    beds = cfg.map { c -> c.central.annotations ? tsvRows(c.central.annotations).collect { r -> resolvePath(r.bed, r._base) } : [] }
    CEN_SUMMARY(features, loci, fm_tables, co_tables, CEN_INTERACTIONS.out.dir, annotations, beds, ingest, exports, config)

    emit:
    summary = CEN_SUMMARY.out.dir
}
