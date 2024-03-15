process execute {
    tag "${data.name}"
    input:
        tuple val(data), path(assembly)
    output:
        tuple val(data.name), path("${assembly}.harvest.combine.scn")
    conda 'bioconda::genometools-genometools'
    cpus 8
    memory 8.GB
    time '6h'
    publishDir 'out_ltr_harvest'
"""
perl ${projectDir}/bin/LTR_HARVEST_parallel/LTR_HARVEST_parallel \
    -seq ${assembly} \
    -threads ${task.cpus} \
    -size 10000000 \
    -time 300
"""
}

workflow LTR_HARVEST {
    take:
        data
    
    main:
        execute(data)       
    
    emit: 
        ltrs = execute.out
}
