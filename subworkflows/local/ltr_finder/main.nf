// version: v1.3

params.chunksize = 5000000
params.timeout = 1500
params.threads = 4
params.params = "-w 2 -C -D 15000 -d 1000 -L 7000 -l 100 -p 20 -M 0.85"

annotator = "LTR_FINDER_parallel_nf"

process chunk_fasta {
    input:
        path(genome)
    output:
        tuple val("${genome.baseName}"), stdout

    """
    fffx chunk ${genome} ${params.chunksize}
    """
}

process ltr_finder {
    input:
        tuple val(name), path(genome_chunk)
    output:
        path("${name}.finder.scn")

    """
    ltr_finder ${params.params} ${genome_chunk} > ${genome}.finder.scn
    """
}

// Next is line 204 in LTR_finder_parallel perl script

workflow LTR_FINDER {
    take:
        genomes
    
    main:
        // Probably not correct either
        chunk_fasta(genomes).map { name, chunk -> tuple(name, ltr_finder(name, chunk)) }
        
}