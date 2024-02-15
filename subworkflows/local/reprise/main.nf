// Has to be compiled manually from source
// TODO: Singularity / docker container

// NOTE: Can be a memory hog, best to leave it disabled for now

process reprise {
    tag "${genome.baseName}"
    input:
        path(genome)
    output:
        path("${genome.baseName}.clustered.reprof")
    cpus 16
    memory '22 GB'
    publishDir 'out_reprise'
    conda 'bioconda::cd-hit bioconda::tantan'

    """
# This step moved to sanitize
#    reformat.sh in=${genome} out=filtered.fa minlength=1000

# Tantan set to mask more (0.005, default, to 0.02)
    tantan -r 0.02 ${genome} > tantan.fa
    nice REPrise  -input tantan.fa -output ${genome.baseName} -pa ${task.cpus} -maxextend 100000 -dist 1
    cd-hit-est -M ${task.memory.toMega()} -T ${task.cpus} -i ${genome.baseName}.reprof -o ${genome.baseName}.clustered.reprof -c 0.8 -p 1
    """
}

workflow REPRISE {
    take:
        genomes
    
    main:
        reprise(genomes)
        
}
