process CLEANUP_TANDEM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/gallvp/edta-components:v0.1"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    if ( "$fasta" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    cleanup_tandem.pl \\
        $args \\
        -f $fasta \\
        > ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's/This is \\(.*\\) built.*/\\1/p')
        trf: \$(trf -v |& sed -n 's/.*Version \\(.*\\)/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    if ( "$fasta" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's/This is \\(.*\\) built.*/\\1/p')
        trf: \$(trf -v |& sed -n 's/.*Version \\(.*\\)/\\1/p')
    END_VERSIONS
    """
}
