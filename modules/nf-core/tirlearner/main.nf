process TIRLEARNER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tir-learner%3A3.0.1--hdfd78af_0':
        'biocontainers/tir-learner%3A3.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa")       , emit: fasta
    tuple val(meta), path("*.gff3")     , emit: gff
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def species         = task.ext.species ?: 'others'
    def maxint          = task.ext.maxint ?: 5000
    """
    TIR-Learner \
        -f ${fasta} \
        -s ${species} \
        -t ${task.cpus} \
        -l ${maxint} \
        -o ${prefix}.TIR \
        ${args}

    mv "${prefix}.TIR/TIR-Learner-Result/TIR-Learner_FinalAnn_filter.fa"    "${prefix}.fa"
    mv "${prefix}.TIR/TIR-Learner-Result/TIR-Learner_FinalAnn_filter.gff3"  "${prefix}.gff3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TIR-Learner: \$(TIR-Learner -v | head -n 1 | sed 's/TIR-Learner //')
    END_VERSIONS
    """

    stub:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.fa"
    touch "${prefix}.gff3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TIR-Learner: \$(TIR-Learner -v | head -n 1 | sed 's/TIR-Learner //')
    END_VERSIONS
    """
}
