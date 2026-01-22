#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    galantelab/precisecaller
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/galantelab/precisecaller
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PRECISECALLER           } from './workflows/precisecaller'
include { PREPARE_GENOME          } from './subworkflows/local/prepare_genome'
include { PREPARE_INTERVALS       } from './subworkflows/local/prepare_intervals'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_precisecaller_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_precisecaller_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_precisecaller_pipeline'
include { getAssayAttribute       } from './subworkflows/local/utils_nfcore_precisecaller_pipeline'
include { paramsHelp              } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT HELP MESSAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info paramsHelp("nextflow run galantelab/precisecaller --input input_file.csv")
    exit 0
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta             = getGenomeAttribute('fasta')
params.fasta_fai         = getGenomeAttribute('fasta_fai')
params.dict              = getGenomeAttribute('dict')
params.bwa               = getGenomeAttribute('bwa')
params.dbsnp             = getGenomeAttribute('dbsnp')
params.dbsnp_tbi         = getGenomeAttribute('dbsnp_tbi')
params.known_indels      = getGenomeAttribute('known_indels')
params.known_indels_tbi  = getGenomeAttribute('known_indels_tbi')
params.known_snps        = getGenomeAttribute('known_snps')
params.known_snps_tbi    = getGenomeAttribute('known_snps_tbi')
params.vep_cache_version = getGenomeAttribute('vep_cache_version')
params.vep_genome        = getGenomeAttribute('vep_genome')
params.vep_species       = getGenomeAttribute('vep_species')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow GALANTELAB_PRECISECALLER {
    take:
    samplesheet

    main:
    versions = Channel.empty()

    //
    // SUBWORKFLOW: Prepare reference genome files
    //
    PREPARE_GENOME(
        params.fasta,
        params.fasta_fai,
        params.dict,
        params.bwa,
        params.dbsnp,
        params.dbsnp_tbi,
        params.known_indels,
        params.known_indels_tbi,
        params.known_snps,
        params.known_snps_tbi,
        params.vep_cache,
        params.vep_cache_version,
        params.vep_genome,
        params.vep_species
    )
    versions = versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Prepare intervals file
    //
    PREPARE_INTERVALS(
        PREPARE_GENOME.out.fasta_fai,
        params.intervals
    )
    versions = versions.mix(PREPARE_INTERVALS.out.versions)

    //
    // Prepare files based on params
    //
    umi_file = params.umi_file ?
        file(params.umi_file, checkIfExists: true) :
        null

    //
    // Prepare filters for variants
    //
    filters_indel_map = GatkFilters.parse(params.filters_indel)
    filters_snp_map   = GatkFilters.parse(params.filters_snp)

    //
    // WORKFLOW: Run pipeline
    //
    PRECISECALLER (
        samplesheet,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fasta_fai,
        PREPARE_GENOME.out.dict,
        PREPARE_GENOME.out.bwa,
        PREPARE_GENOME.out.known_indels,
        PREPARE_GENOME.out.known_indels_tbi,
        PREPARE_GENOME.out.known_snps,
        PREPARE_GENOME.out.known_snps_tbi,
        PREPARE_GENOME.out.vep_cache,
        PREPARE_INTERVALS.out.intervals,
        PREPARE_INTERVALS.out.intervals_gz_tbi,
        umi_file,
        filters_indel_map,
        filters_snp_map,
        versions
    )

    emit:
    multiqc_report = PRECISECALLER.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    GALANTELAB_PRECISECALLER (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        GALANTELAB_PRECISECALLER.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
