include { ANNOSINE                                  } from '../../modules/edta-components/annosine'
include { TESORTER as ANNOSINE_TESORTER             } from '../../modules/edta-components/tesorter/main'
include { CLEANUP_MISCLAS                           } from '../../modules/local/cleanup_misclas/main'
include { CLEANUP_TANDEM as ANNOSINE_CLEANUP_TANDEM } from '../../modules/local/cleanup_tandem/main'

workflow FASTA_ANNOSINE {
    
    take:
    ch_short_ids_fasta              // Channel: [ meta, fasta ]; meta -> [ val(id) ]

    main:
    // Versions
    ch_versions                     = Channel.empty()

    // MODULE: ANNOSINE
    ANNOSINE(
        ch_short_ids_fasta,
        3 // mode
    )

    ch_annosine_fa                  = ANNOSINE.out.fa
    ch_versions                     = ch_versions.mix(ANNOSINE.out.versions.first())

    // MODULE: TESORTER as ANNOSINE_TESORTER
    ANNOSINE_TESORTER(
        ch_annosine_fa,
        [] // db_hmm
    )

    ch_tesorter_tsv                 = ANNOSINE_TESORTER.out.cls_tsv
    ch_versions                     = ch_versions.mix(ANNOSINE_TESORTER.out.versions.first())

    // MODULE: CLEANUP_MISCLAS
    ch_cleanup_inputs               = ch_tesorter_tsv
                                    | join(ch_annosine_fa)

    CLEANUP_MISCLAS(
        ch_cleanup_inputs.map { meta, tsv, fa -> [ meta, tsv ] },
        ch_cleanup_inputs.map { meta, tsv, fa -> fa }
    )

    ch_cleanup_misclas_cln          = CLEANUP_MISCLAS.out.cln
    ch_versions                     = ch_versions.mix(CLEANUP_MISCLAS.out.versions.first())

    // MODULE: CLEANUP_TANDEM as ANNOSINE_CLEANUP_TANDEM
    ANNOSINE_CLEANUP_TANDEM ( ch_cleanup_misclas_cln )

    ch_cleanup_tandem_fasta         = ANNOSINE_CLEANUP_TANDEM.out.fasta
    ch_versions                     = ch_versions.mix(ANNOSINE_CLEANUP_TANDEM.out.versions.first())

    emit:
    fasta                           = ch_cleanup_tandem_fasta   // [ meta, fasta ]
    versions                        = ch_versions               // Channel: [ versions.yml ]
}