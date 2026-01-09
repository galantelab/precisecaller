//
// Prepare reference genome files
//

include { BWA_INDEX                                             } from '../../../modules/nf-core/bwa/index/main'
include { GATK4_CREATESEQUENCEDICTIONARY                        } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GUNZIP                          as GUNZIP_FASTA       } from '../../../modules/nf-core/gunzip/main'
include { BCFTOOLS_SORT                   as SORT_INDEX_DBSNP   } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT                   as SORT_INDEX_INDELS  } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_SORT                   as SORT_INDEX_SNPS    } from '../../../modules/nf-core/bcftools/sort/main'
include { SAMTOOLS_FAIDX                                        } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {
    take:
    fasta                // [mandatory] value: string - filename
    fasta_fai            // [optional]  value: string - filename
    dict                 // [optional]  value: string - filename
    bwa                  // [optional]  value: string - dirname
    dbsnp                // [mandatory] value: string - filename
    dbsnp_tbi            // [optional]  value: string - filename
    known_indels         // [mandatory] value: string - filename
    known_indels_tbi     // [optional]  value: string - filename
    known_snps           // [mandatory] value: string - filename
    known_snps_tbi       // [optional]  value: string - filename

    main:
    ch_versions = Channel.empty()

    fasta_file        = file(fasta, checkIfExists: true)
    dbsnp_file        = file(dbsnp, checkIfExists: true)
    known_indels_file = file(known_indels, checkIfExists: true)
    known_snps_file   = file(known_snps, checkIfExists: true)

    ch_fasta = Channel.empty()
    if (fasta.endsWith('.gz')) {
        GUNZIP_FASTA([[id:"${fasta_file.baseName}"], fasta_file])
        ch_fasta     = GUNZIP_FASTA.out.gunzip.first()
        ch_versions  = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value([[id:"${fasta_file.baseName}"], fasta_file])
    }

    ch_fasta_fai = Channel.empty()
    if (fasta_fai) {
        ch_fasta_fai = Channel.value([[id:"fai"], file(fasta_fai, checkIfExists: true)])
    } else {
        SAMTOOLS_FAIDX(ch_fasta, [[ id:'no_fai' ], []], false)
        ch_fasta_fai = SAMTOOLS_FAIDX.out.fai.first()
        ch_versions  = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    ch_dict = Channel.empty()
    if (dict) {
        ch_dict = Channel.value([[id:"dict"], file(dict, checkIfExists: true)])
    } else {
        GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
        ch_dict     = GATK4_CREATESEQUENCEDICTIONARY.out.dict.first()
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

    ch_dbsnp     = Channel.value([[id:"${dbsnp_file.baseName}"], dbsnp_file])
    ch_dbsnp_tbi = Channel.empty()

    if (dbsnp_tbi) {
        dbsnp_tbi_file = file(dbsnp_tbi, checkIfExists: true)
        ch_dbsnp_tbi   = Channel.value([[id:"${dbsnp_file.baseName}"], dbsnp_tbi_file])
    } else {
        SORT_INDEX_DBSNP(ch_dbsnp)
        ch_dbsnp     = SORT_INDEX_DBSNP.out.vcf.first()
        ch_dbsnp_tbi = SORT_INDEX_DBSNP.out.tbi.first()
        ch_versions  = ch_versions.mix(SORT_INDEX_DBSNP.out.versions)
    }

    ch_known_indels     = Channel.value([[id:"${known_indels_file.baseName}"], known_indels_file])
    ch_known_indels_tbi = Channel.empty()

    if (known_indels_tbi) {
        known_indels_tbi_file = file(known_indels_tbi, checkIfExists: true)
        ch_known_indels_tbi   = Channel.value([[id:"${known_indels_file.baseName}"], known_indels_tbi_file])
    } else {
        SORT_INDEX_INDELS(ch_known_indels)
        ch_known_indels     = SORT_INDEX_INDELS.out.vcf.first()
        ch_known_indels_tbi = SORT_INDEX_INDELS.out.tbi.first()
        ch_versions         = ch_versions.mix(SORT_INDEX_INDELS.out.versions)
    }

    ch_known_snps     = Channel.value([[id:"${known_snps_file.baseName}"], known_snps_file])
    ch_known_snps_tbi = Channel.empty()

    if (known_snps_tbi) {
        known_snps_tbi_file = file(known_snps_tbi, checkIfExists: true)
        ch_known_snps_tbi   = Channel.value([[id:"${known_snps_file.baseName}"], known_snps_tbi_file])
    } else {
        SORT_INDEX_SNPS(ch_known_snps)
        ch_known_snps     = SORT_INDEX_SNPS.out.vcf.first()
        ch_known_snps_tbi = SORT_INDEX_SNPS.out.tbi.first()
        ch_versions       = ch_versions.mix(SORT_INDEX_SNPS.out.versions)
    }

    emit:
    fasta            = ch_fasta             // channel: tuple(meta, fasta)
    fasta_fai        = ch_fasta_fai         // channel: tuple(meta, fasta.fai)
    bwa              = ch_bwa               // channel: tuple(meta, bwa)
    dbsnp            = ch_dbsnp             // channel: tuple(meta, dbsnp)
    dbsnp_tbi        = ch_dbsnp_tbi         // channel: tuple(meta, dbsnp_tbi)
    known_indels     = ch_known_indels      // channel: tuple(meta, known_indels)
    known_indels_tbi = ch_known_indels_tbi  // channel: tuple(meta, known_indels_tbi)
    known_snps       = ch_known_snps        // channel: tuple(meta, known_snps)
    known_snps_tbi   = ch_known_snps_tbi    // channel: tuple(meta, known_snps_tbi)
    versions         = ch_versions          // channel: path(versions.yml)
}
