// Has to be compiled manually from source
// TODO: Singularity / docker container

process reprise {
    tag "${genome.baseName}"
    input:
        path(genome)
    output:
        path("${genome.baseName}.reprof")
    cpus 16
    publishDir 'out_reprise'

    """
    REPrise -input ${genome} -output ${genome.baseName} -pa ${task.cpus}
    """
}

workflow REPRISE {
    take:
        genomes
    
    main:
        reprise(genomes)
        
}