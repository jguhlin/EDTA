process POSTPROCESS_LTRRETRIEVER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/gallvp/edta-components:v0.1"

    input:
    tuple val(meta), path(fasta)
    path(pass_list)
    path(defalse)
    path(pass_list_gff)

    output:
    // TODO: Not sure about other outputs yet
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def MDUST_VERSION = '2006.10.17' // mdust CLI does not report version
    """
    bin_path=\$(dirname \$(which call_seq_by_list.pl))
    
    ln -s \\
        $pass_list \\
        "${fasta}.pass.list"

    ln -s \\
        $defalse \\
        "${fasta}.defalse"

    ln -s \\
        $pass_list_gff \\
        "${fasta}.pass.list.gff3"

    postprocess_ltrretriever.pl \\
        $fasta \\
        $task.cpus \\
        \$bin_path

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's/This is \\(.*\\) built.*/\\1/p')
        mdust: $MDUST_VERSION
        TEsorter: \$(TEsorter -v | tr -d 'TEsorter ')
        trf: \$(trf -v |& sed -n 's/.*Version \\(.*\\)/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def MDUST_VERSION = '2006.10.17' // mdust CLI does not report version
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's/This is \\(.*\\) built.*/\\1/p')
        mdust: $MDUST_VERSION
        TEsorter: \$(TEsorter -v | tr -d 'TEsorter ')
        trf: \$(trf -v |& sed -n 's/.*Version \\(.*\\)/\\1/p')
    END_VERSIONS
    """
}
