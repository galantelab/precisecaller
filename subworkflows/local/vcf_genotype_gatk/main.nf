//
// Joint genotyping of multiple GVCFs using GATK
//
// This subworkflow supports two alternative strategies for combining per-sample GVCFs:
//
// 1) GenomicsDBImport  (recommended for WGS / large cohorts)
//    - Scales efficiently for high coverage and large sample counts
//    - Uses a workspace database rather than producing an intermediate combined GVCF
//
// 2) CombineGVCFs  (recommended for WXS / panels / smaller cohorts)
//    - Produces a merged multi-sample GVCF file directly
//    - Simpler, but less scalable for very large datasets
//
// The strategy is controlled by `use_genomicsdb`.
// Upstream logic (e.g. assay defaults) should decide the appropriate mode
// The final output is always a joint-called VCF + TBI from GenotypeGVCFs
//

include { GATK4_GENOMICSDBIMPORT as GENOMICSDBIMPORT } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_COMBINEGVCFS     as COMBINEGVCFS     } from '../../../modules/nf-core/gatk4/combinegvcfs/main'
include { GATK4_GENOTYPEGVCFS    as GENOTYPEGVCFS    } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'

workflow VCF_GENOTYPE_GATK {
    take:
    vcf                 // [mandatory] channel: tuple(meta, vcf)
    tbi                 // [mandatory] channel: tuple(meta, vcf.tbi)
    intervals           // [mandatory] channel: path(intervals)
    intervals_gz_tbi    // [mandatory] channel: tuple(intervals.gz, intervals_gz_tbi)
    fasta               // [mandatory] channel: tuple(meta, fa)
    fasta_fai           // [mandatory] channel: tuple(meta, fai)
    dict                // [mandatory] channel: tuple(meta, dict)
    dbsnp               // [mandatory] channel: tuple(meta, dbsnp)
    dbsnp_tbi           // [mandatory] channel: tuple(meta, dbsnp_tbi)
    use_genomicsdb      // [optional]  value: boolean

    main:
    versions       = Channel.empty()
    genotype_input = Channel.empty()

    // Group all per-sample GVCFs into a single aggregated tuple.
    // The meta information is replaced by a synthetic meta object
    // identifying that the dataset now represents "all samples"

    grouped_vcf = vcf
        .map { meta, vcf ->
            tuple([id:"all_samples"], vcf)
        }
        .groupTuple()

    grouped_tbi = tbi
        .map { meta, tbi ->
            tuple([id:"all_samples"], tbi)
        }
        .groupTuple()

    // Used primarily for WGS datasets or large cohorts where memory
    // and runtime can become limiting
    if (use_genomicsdb) {
        // These flags configure GenomicsDBImport to operate using BED
        // intervals instead of Picard-style interval lists and without
        // requiring workspace update or input maps.
        // If future use-cases change, these can be parameterized
        interval_value   = false
        run_intlist      = false
        run_updatewspace = false
        input_map        = false

        GENOMICSDBIMPORT(
            grouped_vcf
               .join(grouped_tbi)
               .combine(intervals)
               .map { meta, vcf, tbi, bed ->
                   tuple(meta, vcf, tbi, bed, interval_value, [])
               },
            run_intlist,
            run_updatewspace,
            input_map
        )

        genotype_input = GENOMICSDBIMPORT.out.genomicsdb
                         .combine(intervals_gz_tbi)
                         .map { meta, db, bed, bed_tbi ->
                             tuple(meta, db, [], bed, bed_tbi)
                         }

        versions = versions.mix(GENOMICSDBIMPORT.out.versions)
    } else {
        // Typically preferred for WXS or targeted sequencing panels and
        // for smaller sample sets
        COMBINEGVCFS(
            grouped_vcf.join(grouped_tbi),
            fasta.map     { it[1] },
            fasta_fai.map { it[1] },
            dict.map      { it[1] }
        )

        genotype_input = COMBINEGVCFS.out.combined_gvcf
                         .join(COMBINEGVCFS.out.combined_tbi)
                         .combine(intervals_gz_tbi)
                         .map { meta, vcf, vcf_tbi, bed, bed_tbi ->
                             tuple(meta, vcf, vcf_tbi, bed, bed_tbi)
                         }

        versions = versions.mix(COMBINEGVCFS.out.versions)
    }

    // Regardless of whether GVCFs were combined via GenomicsDBImport
    // or CombineGVCFs, GenotypeGVCFs performs the actual cohort
    // genotyping and produces the final multi-sample VCF output
    GENOTYPEGVCFS(
        genotype_input,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi
    )

    versions = versions.mix(GENOTYPEGVCFS.out.versions)

    emit:
    vcf      = GENOTYPEGVCFS.out.vcf  // channel: tuple(meta, vcf)
    tbi      = GENOTYPEGVCFS.out.tbi  // channel: tuple(meta, vcf.tbi)
    versions = versions               // channel: path(versions.yml)
}
