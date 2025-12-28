//
// Prepare reference genome files
//

include { BWA_INDEX                      } from '../../../modules/nf-core/bwa/index/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GUNZIP as GUNZIP_FASTA         } from '../../../modules/nf-core/gunzip'
include { SAMTOOLS_FAIDX                 } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {
    take:
    fasta      // [mandatory] path(fasta)
    fasta_fai  // [optional]  path(fasta.fai)
    dict       // [optional]  path(dict)
    bwa        // [optional]  path(bwa)

    main:
    ch_versions = Channel.empty()
    fasta_file  = file(fasta, checkIfExists: true)

    ch_fasta = Channel.empty()
    if (fasta.endsWith('.gz')) {
        GUNZIP_FASTA([[id:"${fasta_file.baseName}"], fasta_file])
        ch_fasta     = GUNZIP_FASTA.out.gunzip
        ch_versions  = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value([[id:"${fasta_file.baseName}"], fasta_file])
    }

    ch_fasta_fai = Channel.empty()
    if (fasta_fai) {
        ch_fasta_fai = Channel.value([[id:"fai"], file(fasta_fai, checkIfExists: true)])
    } else {
        SAMTOOLS_FAIDX(ch_fasta, [[ id:'no_fai' ], []], false)
        ch_fasta_fai = SAMTOOLS_FAIDX.out.fai
        ch_versions  = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    ch_dict = Channel.empty()
    if (dict) {
        ch_dict = Channel.value([[id:"dict"], file(dict, checkIfExists: true)])
    } else {
        GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
        ch_dict     = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }

    ch_bwa = Channel.empty()
    if (bwa) {
        ch_bwa = Channel.value([[id:"bwa"], file(bwa, checkIfExists: true)])
    } else {
        BWA_INDEX(ch_fasta)
        ch_bwa       = BWA_INDEX.out.index
        ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
    }

    emit:
    fasta     = ch_fasta      // channel: tuple(meta, fasta)
    fasta_fai = ch_fasta_fai  // channel: tuple(meta, fasta.fai)
    bwa       = ch_bwa        // channel: tuple(meta, bwa)
    versions  = ch_versions   // channel: path(versions.yml)
}
