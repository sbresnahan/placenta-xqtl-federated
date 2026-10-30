/*
 * --mode site : genotype QC -> ancestry strata -> per-stratum genotype QC -> normalization ->
 * HCP grid and k selection -> tensorQTL (permutation, conditionally independent, nominal) ->
 * combined-modality arm -> LD -> export (files A-F, I, J) with the mandatory leak scan.
 *
 * Phenotype inputs come from params.phenotypes_dir: <modality>.bed.gz per modality (gene-TSS
 * coordinates; grouped modalities use GENE__feature IDs), expression.counts.bed.gz (raw counts for
 * TMM-VST), qc_metrics.tsv and cell_proportions.tsv (per RNA-seq sample).
 */

include { GENO_QC ; GENO_MERGE ; ANCESTRY ; STRATUM_GENO ; LD_MATRICES } from '../modules/local/genotypes.nf'
include { NORMALIZE ; NORMALIZE as NORMALIZE_ISOFORM ; COMBINE_MODALITIES ; HCP_EVAL ; SELECT_K ;
          CIS_PERMUTATION ; CIS_INDEPENDENT ; MERGE_INDEPENDENT ; CIS_NOMINAL } from '../modules/local/mapping.nf'
include { VERIFY_RESTRICTED ; EXPORT_STRATUM ; EXPORT_FINALIZE } from '../modules/local/export.nf'

// HCP grid: cohort.yml mapping.hcp_grid, else 0..60 by 10 (n < 350) or 0..100 by 10 (atlas rule).
// k values that leave fewer than 10 residual degrees of freedom are dropped.
def hcpGrid(c, samples) {
    def n = samples.readLines().size() - 1
    def grid = c.mapping.hcp_grid ?: (n < 350 ? (0..60).step(10) : (0..100).step(10))
    return grid.collect { k -> k as int }.findAll { k -> k == 0 || k < n - 20 }
}

def nSignificant(summary) {
    def lines = summary.readLines()
    def hdr = lines[0].tokenize('\t')
    return lines[1].tokenize('\t')[hdr.indexOf('n_significant')] as int
}

def phenoFile(c, name) {
    return file("${c.phenotypes_dir}/${name}")
}

def isIsoform(m) {
    return m in ['isoform_expression', 'isoform_usage']
}

workflow SITE {
    take:
    cfg     // value channel: resolved cohort parameters (map)

    main:
    mcfg = cfg.map { c -> c.mapping + [maf: c.qc.maf_min, exclude_cell_types: c.gxc.exclude_cell_types] }
    metadata = cfg.map { c -> file(c.metadata, checkIfExists: true) }
    qc_metrics = cfg.map { c -> phenoFile(c, 'qc_metrics.tsv') }
    cell_props = cfg.map { c -> phenoFile(c, 'cell_proportions.tsv').exists() ? phenoFile(c, 'cell_proportions.tsv') : [] }

    // ---------------- genotypes ----------------
    vcfs = cfg.flatMap { c ->
        (1..22).collect { n -> tuple("chr${n}".toString(), file("${c.genotypes.imputed_dir}/${c.genotypes.imputed_pattern.replace('{chr}', n.toString())}")) }
            .findAll { t -> t[1].exists() }
    }
        .ifEmpty { error("No imputed VCFs found: check genotypes.imputed_dir and genotypes.imputed_pattern in cohort.yml") }
    GENO_QC(vcfs, metadata, cfg.map { c -> c.qc.rsq_min })
    chroms = vcfs.map { t -> t[0] }.toSortedList { a, b -> (a - 'chr') as int <=> (b - 'chr') as int }
    GENO_MERGE(GENO_QC.out.files.collect(), chroms)

    ANCESTRY(GENO_MERGE.out.pfile,
        cfg.map { c -> files("${c.refs_dir}/ancestry/1kGP_hg38_pca_panel.*") },
        cfg.map { c -> file("${c.refs_dir}/ancestry/1kg_superpop.tsv", checkIfExists: true) },
        metadata, cfg.map { c -> c.ancestry })
    strata = ANCESTRY.out.samples.flatten()
        .map { f -> tuple(f.name.replaceAll(/^stratum_/, '').replaceAll(/\.samples\.tsv$/, ''), f) }

    STRATUM_GENO(strata.combine(GENO_MERGE.out.pfile).combine(GENO_MERGE.out.info), cfg.map { c -> c.qc })
    geno = STRATUM_GENO.out.geno
    LD_MATRICES(geno.join(STRATUM_GENO.out.variant_qc), cfg.map { c -> c.ld })

    // ---------------- phenotypes: normalization within stratum ----------------
    mods = cfg.flatMap { c ->
        c.mapping.modalities.collect { m ->
            def bed = phenoFile(c, "${m}.bed.gz")
            if (!bed.exists()) {
                error("Modality '${m}' is listed in mapping.modalities but ${bed} does not exist")
            }
            def counts = (m == 'expression' && phenoFile(c, 'expression.counts.bed.gz').exists()) ? phenoFile(c, 'expression.counts.bed.gz') : []
            tuple(m, bed, counts)
        }
    }
    norm_in = STRATUM_GENO.out.samples.combine(mods)
        .map { anc, smp, m, bed, counts -> tuple(anc, m, bed, counts, smp) }
    NORMALIZE(norm_in.filter { t -> !isIsoform(t[1]) }.map { anc, m, bed, counts, smp -> tuple(anc, m, bed, counts, [], smp) })
    // isoform modalities drop the samples flagged as expression connectivity outliers
    expr_outliers = NORMALIZE.out.outliers.filter { t -> t[1] == 'expression' }.map { anc, _m, o -> tuple(anc, o) }
    NORMALIZE_ISOFORM(norm_in.filter { t -> isIsoform(t[1]) }.combine(expr_outliers, by: 0)
        .map { anc, m, bed, counts, smp, o -> tuple(anc, m, bed, counts, o, smp) })
    norm_bed = NORMALIZE.out.bed.mix(NORMALIZE_ISOFORM.out.bed)
    norm_staging = NORMALIZE.out.staging.mix(NORMALIZE_ISOFORM.out.staging)

    // combined arm: all modalities of a stratum, one group per gene
    COMBINE_MODALITIES(norm_bed.map { anc, m, bed, _g -> tuple(anc, m, bed) }.groupTuple()
        .filter { t -> t[1].size() > 1 })
    units = norm_bed.mix(COMBINE_MODALITIES.out.bed.map { anc, m, bed, g -> tuple(anc, m, bed, [g]) })

    // ---------------- HCP grid and k* ----------------
    hcp_in = units.combine(STRATUM_GENO.out.samples, by: 0).combine(STRATUM_GENO.out.pcs, by: 0)
        .combine(geno, by: 0).combine(cfg)
        .flatMap { anc, m, bed, g, smp, pcs, pgen, pvar, psam, c ->
            hcpGrid(c, smp).collect { k -> tuple(anc, m, k, bed, g, smp, pcs, pgen, pvar, psam) }
        }
    HCP_EVAL(hcp_in, metadata, qc_metrics, cell_props, mcfg)

    evals = HCP_EVAL.out.eval.map { anc, m, _k, ev, _cov -> tuple(anc, m, ev) }.groupTuple(by: [0, 1])
    sel_in = units.map { anc, m, _bed, g -> tuple(anc, m, g ? true : false) }
        .join(evals, by: [0, 1])
        .join(norm_staging.map { anc, m, _f, ns -> tuple(anc, m, ns) }, by: [0, 1], remainder: true)
        .filter { t -> t[2] != null }
        .map { anc, m, grouped, evs, ns ->
            def arms = (grouped && m != 'combined') ? [m, "${m}_ungrouped".toString()] : [m]
            tuple(anc, m, arms, grouped, ns ?: [], evs)
        }
    SELECT_K(sel_in, mcfg)
    kstar = SELECT_K.out.k.map { anc, m, kf -> tuple(anc, m, kf.text.trim() as int) }
    cov_sel = HCP_EVAL.out.eval.map { anc, m, k, _ev, cov -> tuple(anc, m, k as int, cov) }
        .join(kstar, by: [0, 1, 2])
        .map { anc, m, _k, cov -> tuple(anc, m, cov) }

    // ---------------- arms: grouped + ungrouped per modality, combined ----------------
    unit_cov = units.join(cov_sel, by: [0, 1])
    arms = unit_cov.flatMap { anc, m, bed, g, cov ->
        def out = [tuple(anc, m, bed, g, cov)]
        if (g && m != 'combined') {
            out << tuple(anc, "${m}_ungrouped".toString(), bed, [], cov)
        }
        out
    }
        .combine(geno, by: 0)
    CIS_PERMUTATION(arms, mcfg)

    ind_in = CIS_PERMUTATION.out.cis.combine(cfg)
        .flatMap { anc, arm, cis, summ, c ->
            def nc = Math.max(1, Math.ceil(nSignificant(summ) / (c.mapping.independent_chunk_size as double)) as int)
            (1..nc).collect { i -> tuple(anc, arm, i, cis) }
        }
        .combine(arms, by: [0, 1])
        .map { anc, arm, i, cis, bed, g, cov, pgen, pvar, psam -> tuple(anc, arm, i, bed, g, cov, pgen, pvar, psam, cis) }
    CIS_INDEPENDENT(ind_in, mcfg)
    MERGE_INDEPENDENT(CIS_INDEPENDENT.out.chunk.groupTuple(by: [0, 1]))

    // nominal pass: every pair, per modality, ungrouped
    CIS_NOMINAL(unit_cov.filter { t -> t[1] != 'combined' }
        .map { anc, m, bed, _g, cov -> tuple(anc, m, bed, cov) }.combine(geno, by: 0), mcfg)

    // ---------------- export ----------------
    VERIFY_RESTRICTED(cfg.map { c -> file(c.restricted, checkIfExists: true) }, cfg.map { c -> c.expected_restricted_version })
    staging = STRATUM_GENO.out.staging.map { anc, a, b, s -> tuple(anc, [a, b, s]) }
        .mix(norm_staging.map { anc, _m, f, ns -> tuple(anc, [f, ns]) })
        .mix(SELECT_K.out.settings)
        .mix(CIS_PERMUTATION.out.cis.map { anc, _arm, cis, summ -> tuple(anc, [cis, summ]) })
        .mix(CIS_PERMUTATION.out.cov.map { anc, _arm, cov -> tuple(anc, [cov]) })
        .mix(MERGE_INDEPENDENT.out.independent.map { anc, _arm, ind -> tuple(anc, [ind]) })
        .mix(CIS_NOMINAL.out.nominal.map { anc, _m, d -> tuple(anc, [d]) })
        .mix(LD_MATRICES.out.ld.map { anc, d -> tuple(anc, [d]) })
        .transpose()
        .groupTuple()
    ex = cfg.combine(VERIFY_RESTRICTED.out.version).map { c, v ->
        c.export + [cohort_id: c.cohort_id, image_sha256: c.image_sha256 ?: 'unknown',
                    restricted_version: v.text.trim(), pipeline_version: workflow.manifest.version]
    }
        .first()
    EXPORT_STRATUM(staging, qc_metrics, cell_props, ex)
    ids = cfg.map { c -> [file(c.metadata), file(c.samplesheet)] }.combine(GENO_MERGE.out.vcf_samples).first()
    EXPORT_FINALIZE(EXPORT_STRATUM.out.dir.collect(), ANCESTRY.out.strata, ids, ex)

    emit:
    export = EXPORT_FINALIZE.out.export
}
