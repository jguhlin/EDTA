// Has to be compiled manually from source
// TODO: Singularity / docker container

process reprise {
    tag "${genome.baseName}"
    input:
        path(genome)
    output:
        path("${genome.baseName}.clustered.reprof")
    cpus 16
    memory '2 GB'
    publishDir 'out_reprise'

    """
    #REPrise -input ${genome} -output ${genome.baseName} -pa ${task.cpus}
    REPrise  -input ${genome} -output ${genome.baseName} -pa ${task.cpus} -maxextend 100000 -dist 1
    cd-hit-est -M ${task.memory.toMega()} -T ${task.cpus} -i ${genome.baseName}.reprof -o ${genome.baseName}.clustered.reprof -c 0.8 -p 1
    """
}

workflow REPRISE {
    take:
        genomes
    
    main:
        reprise(genomes)
        
}
