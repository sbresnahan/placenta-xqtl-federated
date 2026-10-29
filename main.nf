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
      --validate_only           check inputs and stop
      -profile                  site/executor profile, e.g. seadragon | slurm,apptainer | local,apptainer
      -resume                   reuse finished tasks

    Resources and scheduler: docs/HPC.md
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

        if (params.validate_only) {
            VALIDATE.out.report.view { r -> "Validation passed: ${r}" }
        }
        else {
            notYet(params.mode)
        }
    }
}
