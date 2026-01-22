//
// Annotates normalized and indexed VCFs using Ensembl VEP
//
// This subworkflow is a thin orchestration layer around the nf-core
// `ENSEMBLVEP_VEP` module, handling input channel alignment and
// exposing annotated VCF outputs.
//
// Notes:
// - This subworkflow assumes that VCFs have already been normalized
//   (e.g. split multiallelic sites and left-aligned INDELs)
// - The VEP cache is expected to be pre-built and compatible with the
//   provided genome, species, and cache version
// - Version information and MultiQC reports are emitted by the module
//   via topic-based channels and must be handled upstream by the
//   parent workflow

include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep/main'

workflow VCF_ANNOTATION_VEP {
    take:
    vcf                  // [mandatory] channel: tuple(meta, vcf)
    tbi                  // [mandatory] channel: tuple(meta, vcf.tbi)
    fasta                // [mandatory] channel: tuple(meta, fa)
    vep_cache            // [mandatory] channel: tuple(meta, vep_cache/)
    vep_cache_version    // [mandatory] value: string
    vep_genome           // [mandatory] value: string
    vep_species          // [mandatory] value: string

    main:
    ENSEMBLVEP_VEP(
        vcf.join(tbi),
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache.map { it[1] },
        fasta,
        []
    )

    emit:
    vcf  = ENSEMBLVEP_VEP.out.vcf  // channel: tuple(meta, vcf)
    tbi  = ENSEMBLVEP_VEP.out.tbi  // channel: tuple(meta, vcf.tbi)
    // multiqc_files from topic channel
    // versions from topic channel
}
