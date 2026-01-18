// Split multiallelic variants on GVCFs

include { BCFTOOLS_NORM as SPLIT_MULTIALLELIC_SITES } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as LEFT_ALIGN_AND_NORMALIZE } from '../../../modules/nf-core/bcftools/norm/main'

workflow VCF_NORM_BCFTOOLS {
    take:
    vcf    // [mandatory] channel: tuple(meta, vcf)
    tbi    // [mandatory] channel: tuple(meta, vcf.tbi)
    fasta  // [mandatory] channel: tuple(meta, fasta)

    main:
    versions = Channel.empty()

    SPLIT_MULTIALLELIC_SITES(
        vcf.join(tbi),
        fasta
    )

    LEFT_ALIGN_AND_NORMALIZE(
        SPLIT_MULTIALLELIC_SITES.out.vcf
            .join(SPLIT_MULTIALLELIC_SITES.out.tbi),
        fasta
    )

    versions = versions.mix(SPLIT_MULTIALLELIC_SITES.out.versions)
    versions = versions.mix(LEFT_ALIGN_AND_NORMALIZE.out.versions)

    emit:
    vcf      = LEFT_ALIGN_AND_NORMALIZE.out.vcf   // channel: tuple(meta, vcf)
    tbi      = LEFT_ALIGN_AND_NORMALIZE.out.tbi   // channel: tuple(meta, vcf.tbi)
    versions = versions                           // channel: path(versions.yml)
}
