params.chunksize = 5000000
params.timeout = 1500
params.threads = 4
params.params = "-w 2 -C -D 15000 -d 1000 -L 7000 -l 100 -p 20 -M 0.85"

annotator = "LTR_FINDER_parallel_nf"
// TODO: Rustize LTR_FINDER_parallel

process ltr_finder {
    tag "${data.name}"
    input:
        tuple val(data), path(assembly)
    output:
        tuple val(data.name), path("${assembly}.finder.combine.scn")
    cpus 4
    memory 3.GB
    time '2h'
    publishDir 'out_ltr_finder'
"""
perl ${projectDir}/bin/LTR_FINDER_parallel/LTR_FINDER_parallel \
    -seq ${assembly} \
    -threads ${task.cpus} \
    -harvest_out \
    -size 10000000 \
    -time 300
"""
}

workflow LTR_FINDER {
    take:
        genomes
    
    main:
        genomes | ltr_finder
    
    emit: 
        ltr_finder.out
}