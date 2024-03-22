include { LTRHARVEST                } from '../../modules/nf-core/ltrharvest/main'
include { LTRFINDER                 } from '../../modules/nf-core/ltrfinder/main'
include { CAT_CAT                   } from '../../modules/nf-core/cat/cat/main'
include { LTRRETRIEVER_LTRRETRIEVER } from '../../modules/nf-core/ltrretriever/ltrretriever/main'
include { CUSTOM_RESTOREGFFIDS      } from '../../modules/pfr/custom/restoregffids/main'
include { POSTPROCESS_LTRRETRIEVER  } from '../../modules/local/postprocess_ltrretriever/main'

workflow FASTA_LTRRETRIEVER {
    
    take:
    ch_short_ids_fasta              // Channel: [ meta, fasta ]
    ch_short_ids_tsv                // Channel: [ meta, tsv ]

    main:
    // Versions
    ch_versions                     = Channel.empty()


    // MODULE: LTRHARVEST
    LTRHARVEST ( ch_short_ids_fasta )

    ch_ltrharvest_scn               = LTRHARVEST.out.scn
    ch_versions                     = ch_versions.mix(LTRHARVEST.out.versions.first())

    // MODULE: LTRFINDER
    LTRFINDER ( ch_short_ids_fasta )

    ch_ltrfinder_scn                = LTRFINDER.out.scn
    ch_versions                     = ch_versions.mix(LTRFINDER.out.versions.first())

    // MODULE: CAT_CAT
    ch_cat_cat_inputs               = ch_ltrharvest_scn
                                    | join(ch_ltrfinder_scn)
                                    | map { meta, harvested, found -> [ meta, [ harvested, found ] ] }

    CAT_CAT ( ch_cat_cat_inputs )

    ch_ltr_candidates               = CAT_CAT.out.file_out
    ch_versions                     = ch_versions.mix(CAT_CAT.out.versions.first())

    // MODULE: LTRRETRIEVER_LTRRETRIEVER
    ch_ltrretriever_inputs          = ch_short_ids_fasta.join(ch_ltr_candidates)

    LTRRETRIEVER_LTRRETRIEVER (
        ch_ltrretriever_inputs.map { meta, fasta, ltr -> [ meta, fasta ] },
        ch_ltrretriever_inputs.map { meta, fasta, ltr -> ltr },
        [],
        [],
        []
    )

    ch_pass_list                    = LTRRETRIEVER_LTRRETRIEVER.out.pass_list
    ch_pass_list_gff                = LTRRETRIEVER_LTRRETRIEVER.out.pass_list_gff
    ch_ltrlib                       = LTRRETRIEVER_LTRRETRIEVER.out.ltrlib
    ch_defalse                      = LTRRETRIEVER_LTRRETRIEVER.out.defalse
    ch_annotation_out               = LTRRETRIEVER_LTRRETRIEVER.out.annotation_out
    ch_annotation_gff               = LTRRETRIEVER_LTRRETRIEVER.out.annotation_gff
    ch_versions                     = ch_versions.mix(LTRRETRIEVER_LTRRETRIEVER.out.versions.first())

    // MODULE: CUSTOM_RESTOREGFFIDS
    ch_restorable_gff_tsv           = ch_annotation_gff.join(ch_short_ids_tsv)

    CUSTOM_RESTOREGFFIDS (
        ch_restorable_gff_tsv.map { meta, gff, tsv -> [ meta, gff ] },
        ch_restorable_gff_tsv.map { meta, gff, tsv -> tsv }
    )

    ch_restored_gff                 = ch_annotation_gff
                                    | join(CUSTOM_RESTOREGFFIDS.out.restored_ids_gff3, by:0, remainder:true)
                                    | map { meta, gff, restored_gff -> [ meta, restored_gff ?: gff ] }

    ch_versions                     = ch_versions.mix(CUSTOM_RESTOREGFFIDS.out.versions.first())

    // MODULE: POSTPROCESS_LTRRETRIEVER
    ch_postprocess_inputs           = ch_short_ids_fasta
                                    | join(ch_pass_list)
                                    | join(ch_defalse)
                                    | join(ch_pass_list_gff)

    POSTPROCESS_LTRRETRIEVER (
        ch_postprocess_inputs.map { meta, fasta, list, defalse, gff -> [ meta, fasta ] },
        ch_postprocess_inputs.map { meta, fasta, list, defalse, gff -> list },
        ch_postprocess_inputs.map { meta, fasta, list, defalse, gff -> defalse },
        ch_postprocess_inputs.map { meta, fasta, list, defalse, gff -> gff }
    )

    ch_versions                     = ch_versions.mix(POSTPROCESS_LTRRETRIEVER.out.versions.first())

    emit:
    pass_list                       = ch_pass_list          // Channel: [ meta, pass.list ]
    annotation_out                  = ch_annotation_out     // Channel: [ meta, out ]
    annotation_gff                  = ch_restored_gff       // Channel: [ meta, gff ]
    ltrlib                          = ch_ltrlib             // Channel: [ meta, fasta ]
    versions                        = ch_versions           // Channel: [ versions.yml ]
}