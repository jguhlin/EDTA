process repeatmodeler {
    tag "${data.name}"
    input:
        tuple val(data), path(assembly)
    output:
        tuple val(data), path(assembly), path("${data.name}-families.fa")
    conda 'bioconda::repeatmodeler=2.0.5-0'
    cpus 48
    time '30d'
    memory 16.GB

shell:
'''
export BLAST_USAGE_REPORT=false
BuildDatabase -name !{data.name} !{assembly}
RepeatModeler -engine ncbi -threads !{task.cpus} -database !{data.name}
'''
}

process process_output {
    tag "${data.name}"
    input:
        tuple val(data), path(assembly), path(families)
    output:
        tuple val(data.name), path("${data.name}.LINE.raw.fa"), path("${data.name}.RM2.fa"), path("${data.name}.RM2.raw.fa")
    conda 'tesorter=1.4.6=pyhdfd78af_0 bioconda::trf'
    cpus 2
    time '1h'
    memory 8.GB
    publishDir 'out_repeatmodeler'

shell:
'''
awk '{print \$1}' !{families} > !{data.name}.RM2.raw.fa
TEsorter !{data.name}.RM2.raw.fa --disable-pass2 -p !{task.cpus}

perl !{projectDir}/util/cleanup_misclas.pl !{data.name}.RM2.raw.fa.rexdb.cls.tsv

TEsorter !{data.name}.RM2.raw.fa.cln --disable-pass2 -p !{task.cpus}
perl -nle 's/>\\S+\\s+/>/; print \$_' !{data.name}.RM2.raw.fa.cln.rexdb.cls.lib > !{data.name}.RM2.raw.fa.cln

perl !{projectDir}/util/cleanup_tandem.pl -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 80 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f !{data.name}.RM2.raw.fa.cln > !{data.name}.RM2.raw.fa.cln2

grep -P 'LINE|SINE' !{data.name}.RM2.raw.fa.cln2 | perl !{projectDir}/util/output_by_list.pl 1 !{data.name}.RM2.raw.fa.cln2 1 - -FA > !{data.name}.LINE.raw.fa
grep -P 'LINE|SINE' !{data.name}.RM2.raw.fa.cln2 | perl !{projectDir}/util/output_by_list.pl 1 !{data.name}.RM2.raw.fa.cln2 1 - -FA -ex > !{data.name}.RM2.fa
'''
}

workflow REPEATMODELER {
    take:
        genomes

    main:
        genomes | repeatmodeler | process_output

    emit:
        process_output.out
}
