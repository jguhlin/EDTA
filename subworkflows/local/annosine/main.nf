process execute {
    tag "${data.name}"
    input:
        tuple val(data), path(assembly)
    output:
        tuple val(data.name), path("${data.name}.SINE.raw.fa")
    conda 'python=3.10 bioconda::annosine2 bioconda::tesorter biopython'
    cpus 8
    memory 16.GB
    time '4h'
    publishDir 'out_annosine'

shell:
'''
touch Seed_SINE.fa

AnnoSINE_v2 -t !{task.cpus} \
  -a 2 \
  --num_alignments 50000 \
  -rpm 0 \
  --copy_number 3 \
  --shift 100 \
  -auto 1 3 !{data.assembly} .

  if [ -s "Seed_SINE.fa" ]; then
    awk '{print \$1}' Seed_SINE.fa > !{data.name}.AnnoSINE.raw.fa
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
    else
        touch !{data.name}.SINE.raw.fa
    fi

'''
}

workflow ANNOSINE {
    take:
        genomes

    main:
        genomes | execute

    emit:
        execute.out
}