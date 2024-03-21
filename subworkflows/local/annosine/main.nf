process annosine {
    tag "${data.name}"
    input:
        tuple val(data), path(assembly)
    output:
        tuple val(data), path(assembly), path("Seed_SINE.fa"), optional: true
    conda 'bioconda::annosine2'
    cpus 8
    memory 12.GB
    time '12h'

shell:
'''
AnnoSINE_v2 -t !{task.cpus} \
  -a 2 \
  --num_alignments 50000 \
  -rpm 0 \
  --copy_number 3 \
  --shift 100 \
  -auto 1 3 !{data.assembly} .
'''
}

process process_output {
    tag "${data.name}"
    input:
        tuple val(data), path(assembly), path(results)
    output:
        tuple val(data.name), path("${data.name}.SINE.raw.fa")
    conda 'python=3.10 biopython bioconda::tesorter bioconda::trf'
    cpus 2
    memory 4.GB
    time '60m'
    publishDir 'out_annosine'

shell:
'''
awk '{print \$1}' !{results} > !{data.name}.AnnoSINE.raw.fa
TEsorter !{data.name}.AnnoSINE.raw.fa --disable-pass2 -p !{task.cpus}
touch !{data.name}.AnnoSINE.raw.fa.rexdb.cls.tsv
perl !{projectDir}/util/cleanup_misclas.pl !{data.name}.AnnoSINE.raw.fa.rexdb.cls.tsv

perl !{projectDir}/util/cleanup_tandem.pl -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 100 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f !{data.name}.AnnoSINE.raw.fa.cln > !{data.name}.SINE.raw.fa
'''
}

workflow ANNOSINE {
    take:
        genomes

    main:
        genomes | annosine | process_output

    emit:
        process_output.out
}