// Has to be compiled manually from source
// TODO: Singularity / docker container

process reprise {
    input:
        path(genome)
    output:
        path("${genome.baseName}.reprof")
    cpus 16

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