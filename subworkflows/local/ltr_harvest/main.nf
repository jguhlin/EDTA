process execute {
    tag "${genome.baseName}"
    input:
        path(genome)
    output:
        path("${genome}.harvest.combine.scn")
    conda 'bioconda::genometools-genometools'
    cpus 8
    memory 8.GB
    time '6h'
    publishDir 'out_ltr_harvest'
"""
perl ${projectDir}/bin/LTR_HARVEST_parallel/LTR_HARVEST_parallel \
    -seq ${genome} \
    -threads ${task.cpus} \
    -size 10000000 \
    -time 300
"""
}

workflow LTR_HARVEST {
    take:
        genomes
    
    main:
        execute(genomes)
    
    emit: 
        ltrs = execute.out
}
