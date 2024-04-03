process POSTPROCESS_REPEATMODELER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/gallvp/edta-components:v0.1"

    input:
    tuple val(meta), path(fasta)

    output:
    // TODO: Not sure about other outputs yet
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bin_path=\$(dirname \$(which cleanup_tandem.pl))

    ln -s \\
        $fasta \\
        $fasta-families.fa

    postprocess_repeatmodeler.pl \\
        $fasta \\
        $task.cpus \\
        \$bin_path

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's/This is \\(.*\\) built.*/\\1/p')
        TEsorter: \$(TEsorter -v | tr -d 'TEsorter ')
        trf: \$(trf -v |& sed -n 's/.*Version \\(.*\\)/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's/This is \\(.*\\) built.*/\\1/p')
        TEsorter: \$(TEsorter -v | tr -d 'TEsorter ')
        trf: \$(trf -v |& sed -n 's/.*Version \\(.*\\)/\\1/p')
    END_VERSIONS
    """
}
