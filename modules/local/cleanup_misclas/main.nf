process CLEANUP_MISCLAS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/gallvp/edta-components:v0.1"

    input:
    tuple val(meta), path(tsv)
    path fasta

    output:
    tuple val(meta), path('*.cln')          , emit: cln
    tuple val(meta), path('*.cln.list')     , emit: cln_list
    tuple val(meta), path('*.dirt.list')    , emit: dirt_list
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cleanup_misclas.pl \\
        $tsv \\
        $fasta
    
    mv \\
        "${fasta}.cln" \\
        "${prefix}.cln"
    
    mv \\
        "${fasta}.cln.list" \\
        "${prefix}.cln.list"
    
    mv \\
        "${fasta}.dirt.list" \\
        "${prefix}.dirt.list"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's/This is \\(.*\\) built.*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.cln"
    touch "${prefix}.cln.list"
    touch "${prefix}.dirt.list"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's/This is \\(.*\\) built.*/\\1/p')
    END_VERSIONS
    """
}
