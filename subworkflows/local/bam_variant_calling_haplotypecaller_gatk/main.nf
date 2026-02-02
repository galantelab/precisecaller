// GATK HaplotypeCaller calls variants (SNPs and INDELs) via local de-novo
// reassembly of haplotypes in regions showing signs of variation

include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'

workflow BAM_VARIANT_CALLING_HAPLOTYPECALLER_GATK {
    take:
    bam         // [mandatory] channel: tuple(meta, bam)
    bai         // [mandatory] channel: tuple(meta, bai)
    intervals   // [mandatory] channel: path(intervals)
    fasta       // [mandatory] channel: tuple(meta, fa)
    fasta_fai   // [mandatory] channel: tuple(meta, fai)
    dict        // [mandatory] channel: tuple(meta, dict)
    dbsnp       // [mandatory] channel: tuple(meta, dbsnp)
    dbsnp_tbi   // [mandatory] channel: tuple(meta, dbsnp_tbi)

    main:
    versions = Channel.empty()

    GATK4_HAPLOTYPECALLER(
        bam.join(bai)
           .combine(intervals)
           .map { meta, bam, bai, bed ->
               tuple(meta, bam, bai, bed, [])
           },
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi
    )

    versions = versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

    emit:
    vcf      = GATK4_HAPLOTYPECALLER.out.vcf  // channel: tuple(meta, vcf)
    tbi      = GATK4_HAPLOTYPECALLER.out.tbi  // channel: tuple(meta, vcf.tbi)
    versions = versions                       // channel: path(versions.yml)
}
