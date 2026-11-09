#!/usr/bin/env nextflow
/*
 * placenta-xqtl-federated
 * Federated multi-ancestry placental xQTL analysis: site-local mapping and summary-level export.
 * https://github.com/sbresnahan/placenta-xqtl-federated   (GPL-3.0)
 *
 * Usage (sites normally call this through run_xqtl.sh, written by the image's bootstrap):
 *   run_xqtl.sh --mode <mode> -params-file cohort.yml -profile <site>
 * Modes: prepare_references | imputation_prep | site | phenotyping_only | export | central | smoke
 * Documentation: docs/ (rendered at https://sbresnahan.github.io/placenta-xqtl-federated/)
 */

include { SMOKE } from './workflows/smoke.nf'
include { SITE } from './workflows/site.nf'
include { EXPORT_ONLY } from './workflows/export.nf'
include { PREPARE_REFERENCES } from './workflows/prepare_references.nf'
include { PHENOTYPING ; PHENOTYPING_ONLY } from './workflows/phenotyping.nf'
include { IMPUTATION_PREP } from './workflows/imputation_prep.nf'
include { CENTRAL } from './workflows/central.nf'

def modes() {
    return ['prepare_references', 'imputation_prep', 'site', 'phenotyping_only', 'export', 'central', 'smoke']
}

def helpMessage() {
    return """
    placenta-xqtl-federated ${workflow.manifest.version}

    run_xqtl.sh --mode <mode> -params-file cohort.yml -profile <site profile> [options]

    Modes
      smoke               validate a new site (scheduler, container, memory, arrays, tensorQTL)
      prepare_references  build STAR / Salmon / IRFinder-S indices from the embedded references
      imputation_prep     QC and split array/WES/WGS genotypes into TOPMed-ready VCFs
      site                full site analysis: phenotyping, mapping, G×C, GxE, LD, export
      phenotyping_only    phenotype BEDs and normalized matrices only (no export)
      export              (re)assemble the export files A–J from a finished site run
      central             MDACC meta-analysis of exports from several sites

    Common options
      -params-file cohort.yml   cohort settings (see docs/CONFIGURATION.md)
      --outdir DIR              output directory (default: results)
      --restricted FILE         xqtl-restricted-<version>.tar.gz (supplied by email; see docs/SETUP.md)
      --site_refs DIR           indices written by --mode prepare_references
      --phenotypes_dir DIR      map ready phenotype BEDs instead of phenotyping FASTQs (--mode site)
      --validate_only           check inputs and stop
      -profile                  site/executor profile, e.g. seadragon | slurm,apptainer | local,apptainer
      -resume                   reuse finished tasks

    Phenotyping: docs/PHENOTYPING.md      Resources and scheduler: docs/HPC.md
    """.stripIndent()
}

process VALIDATE {
    label 'process_single'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_mode

    input:
    path params_json
    path schema
    val mode

    output:
    path 'resolved_params.json', emit: resolved
    path 'validation_report.tsv', emit: report

    script:
    def nofc = workflow.profile.tokenize(',').contains('test') ? '--no-file-check' : ''
    """
    xqtl_validate.py --params ${params_json} --schema ${schema} --mode ${mode} \\
        --out-json resolved_params.json --out-report validation_report.tsv ${nofc}
    """
}

def banner() {
    return "placenta-xqtl-federated ${workflow.manifest.version} | mode=${params.mode} | cohort=${params.cohort_id ?: 'NA'} | profile=${workflow.profile} | image=${params.sif ?: params.image}"
}

def notYet(String m) {
    error("--mode ${m} is not available in this pipeline version (${workflow.manifest.version}). See CHANGELOG.md.")
}

workflow {
    main:
    if (params.help) {
        log.info(helpMessage())
    }
    else if (!params.mode || !(params.mode in modes())) {
        log.info(helpMessage())
        error("--mode must be one of: ${modes().join(', ')} (got: ${params.mode})")
    }
    else if (params.mode == 'smoke') {
        log.info(banner())
        SMOKE()
    }
    else {
        log.info(banner())
        // Every other mode starts by validating cohort.yml and the input tables.
        def pjson = file("${workflow.workDir}/params_${workflow.sessionId}.json")
        pjson.parent.mkdirs()
        pjson.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(new LinkedHashMap(params)))
        VALIDATE(channel.value(pjson), channel.value(file("${projectDir}/assets/schema/cohort.schema.json")), params.mode)

        // resolved parameters (schema defaults filled in) drive the workflows
        cfg = VALIDATE.out.resolved.map { f -> new groovy.json.JsonSlurperClassic().parseText(f.text) }

        if (params.validate_only) {
            VALIDATE.out.report.view { r -> "Validation passed: ${r}" }
        }
        else if (params.mode == 'prepare_references') {
            PREPARE_REFERENCES(cfg)
            PREPARE_REFERENCES.out.manifest.view { _m -> "Site indices written to ${params.site_refs} (pass --site_refs ${params.site_refs} to --mode site)" }
        }
        else if (params.mode == 'imputation_prep') {
            IMPUTATION_PREP(cfg)
            IMPUTATION_PREP.out.checklist.view { _f -> "TOPMed-ready VCFs and TOPMED_CHECKLIST.md written to ${params.outdir}/imputation_prep (upload manually)" }
        }
        else if (params.mode == 'phenotyping_only') {
            PHENOTYPING_ONLY(cfg)
            PHENOTYPING_ONLY.out.normalized.view { _n -> "Phenotype BEDs: ${params.outdir}/phenotypes; normalized matrices: ${params.outdir}/normalized" }
        }
        else if (params.mode == 'site') {
            // a previous export is removed first, so a failed leak scan never leaves a stale export behind
            file("${params.outdir}/export").deleteDir()
            if (params.phenotypes_dir) {
                SITE(cfg, cfg.map { c -> file(c.phenotypes_dir, checkIfExists: true) })
            }
            else {
                PHENOTYPING(cfg)
                SITE(cfg, PHENOTYPING.out.dir)
            }
            SITE.out.export.view { _d -> "Export written and leak-scanned: ${params.outdir}/export (share only this folder)" }
        }
        else if (params.mode == 'export') {
            file("${params.outdir}/export").deleteDir()
            EXPORT_ONLY(cfg)
            EXPORT_ONLY.out.export.view { _d -> "Export rebuilt and leak-scanned: ${params.outdir}/export (share only this folder)" }
        }
        else if (params.mode == 'central') {
            CENTRAL(cfg, VALIDATE.out.resolved)
            CENTRAL.out.summary.view { _d -> "Central meta-analysis written to ${params.outdir}/central (summary: central/summary/central_summary.tsv)" }
        }
        else {
            notYet(params.mode)
        }
    }
}
