//
// Split multi-sample genotyped GVCF into SNPs & INDELs and apply hard-filter
//
// This subworkflow performs hard-filtering of variants using GATK
//
// The workflow:
//   1) Splits a multi-sample genotyped GVCF into SNPs and INDELs
//      using GATK SelectVariants
//   2) Applies user-defined hard filters to each variant type
//      using GATK VariantFiltration
//
// SNP and INDEL filters are provided independently as Maps:
//   [filter_name: filter_expression]
//
// This design allows:
//   - Independent filtering strategies for SNPs and INDELs
//   - A single unified workflow for both variant types
//   - Clear propagation of variant type information via metadata
//

include { GATK4_SELECTVARIANTS    as GATK4_SELECTVARIANTS_BY_TYPE } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_VARIANTFILTRATION                                 } from '../../../modules/nf-core/gatk4/variantfiltration/main'

workflow VCF_VARIANTFILTRATION_GATK {
    take:
    vcf             // [mandatory] channel: tuple(meta, vcf)
    tbi             // [mandatory] channel: tuple(meta, vcf.tbi)
    intervals       // [mandatory] channel: path(intervals)
    fasta           // [mandatory] channel: tuple(meta, fasta)
    fasta_fai       // [mandatory] channel: tuple(meta, fai)
    dict            // [mandatory] channel: tuple(meta, dict)
    filters_indel   // [optional]  value: map - [filter_name: filter_expr]
    filters_snp     // [optional]  value: map - [filter_name: filter_expr]

    main:
    versions = Channel.empty()

    variant_type = Channel.of(
        ['indel', filters_indel],
        ['snp',   filters_snp]
    )

    // For each input VCF:
    //   - Join with its index and genomic intervals
    //   - Duplicate the stream for SNPs and INDELs
    //   - Annotate metadata with the variant type
    GATK4_SELECTVARIANTS_BY_TYPE(
        vcf
            .join(tbi)
            .combine(intervals)
            .combine(variant_type)
            .map { meta, vcf, tbi, bed, type, filters ->
                def id = "${meta.id}.${type}"
                tuple(meta + [id: id, type: type], vcf, tbi, bed)
            }
    )

    // Each variant type (SNP / INDEL) receives its own filter Map
    // Filters are passed as structured data and translated into
    // GATK-compatible command-line arguments inside the module
    GATK4_VARIANTFILTRATION(
        GATK4_SELECTVARIANTS_BY_TYPE.out.vcf
            .join(GATK4_SELECTVARIANTS_BY_TYPE.out.tbi)
            .map { meta, vcf, tbi ->
                tuple(meta.type, meta, vcf, tbi)
            }
            .join(variant_type)
            .map { type, meta, vcf, tbi, filters ->
                tuple(meta, vcf, tbi, filters)
            },
        fasta,
        fasta_fai,
        dict,
        [[id:"none"], []]
    )

    versions = versions.mix(GATK4_SELECTVARIANTS_BY_TYPE.out.versions)
    versions = versions.mix(GATK4_VARIANTFILTRATION.out.versions)

    emit:
    vcf       = GATK4_VARIANTFILTRATION.out.vcf   // channel: tuple(meta, vcf)
    tbi       = GATK4_VARIANTFILTRATION.out.tbi   // channel: tuple(meta, vcf.tbi)
    versions  = versions                          // channel: path(versions.yml)
}
