include { REPEATMODELER_BUILDDATABASE       } from '../../modules/pfr/repeatmodeler/builddatabase/main'
include { REPEATMODELER_REPEATMODELER       } from '../../modules/pfr/repeatmodeler/repeatmodeler/main'
include { POSTPROCESS_REPEATMODELER         } from '../../modules/local/postprocess_repeatmodeler/main'

workflow FASTA_REPEATMODELER {
    
    take:
    ch_short_ids_fasta              // Channel: [ meta, fasta ]; meta -> [ val(id) ]

    main:
    // Versions
    ch_versions                     = Channel.empty()


    // MODULE: REPEATMODELER_BUILDDATABASE
    REPEATMODELER_BUILDDATABASE ( ch_short_ids_fasta )

    ch_db                           = REPEATMODELER_BUILDDATABASE.out.db
    ch_versions                     = ch_versions.mix(REPEATMODELER_BUILDDATABASE.out.versions.first())

    // MODULE: REPEATMODELER_REPEATMODELER
    REPEATMODELER_REPEATMODELER ( ch_db )

    ch_rm_fasta                     = REPEATMODELER_REPEATMODELER.out.fasta
    ch_versions                     = ch_versions.mix(REPEATMODELER_REPEATMODELER.out.versions.first())

    // MODULE: POSTPROCESS_REPEATMODELER
    POSTPROCESS_REPEATMODELER ( ch_rm_fasta )

    ch_versions                     = ch_versions.mix(POSTPROCESS_REPEATMODELER.out.versions.first())

    emit:
    versions                        = ch_versions           // Channel: [ versions.yml ]
}