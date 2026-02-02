//
// Separate hard-filtered, multi-sample SNP/INDEL GVCFs by sample
//
// This subworkflow takes a joint, hard-filtered multi-sample GVCF (either SNPs or INDELs)
// and splits it back into per-sample GVCFs using GATK SelectVariants
//
// It assumes that the input VCF was produced by a joint genotyping step (e.g. GenotypeGVCFs)
// and that the original sample identifiers are preserved in `meta.samples`
//
// Each output VCF contains variants for a single sample only, while preserving the
// variant type (SNP or INDEL) in the metadata.
//

include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_BY_SAMPLE } from '../../../modules/nf-core/gatk4/selectvariants/main'

workflow VCF_SELECTVARIANTS_BY_SAMPLE_GATK {
    take:
    vcf        // [mandatory] channel: tuple(meta, vcf)
    tbi        // [mandatory] channel: tuple(meta, vcf.tbi)
    intervals  // [mandatory] channel: path(intervals)

    main:
    versions = Channel.empty()

    //
    // For each input joint GVCF:
    //   - iterate over all sample IDs stored in meta.samples
    //   - replicate the same VCF as input
    //   - update metadata to reflect the current sample
    //
    // This expands a single multi-sample VCF into one input tuple per sample,
    // allowing GATK SelectVariants to be executed independently for each sample
    //
    GATK4_SELECTVARIANTS_BY_SAMPLE(
        vcf
            .join(tbi)
            .combine(intervals)
            .flatMap { meta, vcf, tbi, bed ->
                meta.samples.collect { sample_id ->
                    def sample_meta = meta + [
                        sample : sample_id,
                        id     : "${sample_id}.${meta.type}"
                    ]
                    tuple(sample_meta, vcf, tbi, bed)
                }
            }
    )

    versions = versions.mix(GATK4_SELECTVARIANTS_BY_SAMPLE.out.versions)

    emit:
    vcf      = GATK4_SELECTVARIANTS_BY_SAMPLE.out.vcf   // channel: tuple(meta, vcf)
    tbi      = GATK4_SELECTVARIANTS_BY_SAMPLE.out.tbi   // channel: tuple(meta, vcf.tbi)
    versions = versions                                 // channel: path(versions.yml)
}
