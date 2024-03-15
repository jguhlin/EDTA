// Output files
// !{name.baseName}.LTR.intact.raw.fa
// !{name.baseName}.LTRlib.fa
// !{name.baseName}.LTR.intact.fa // seems to be optional
// !{name.baseName}.LTR.intact.gff3 // seems to be optional
// !{name.baseName}.LTR.intact.raw.fa.anno.list

params.u                = 1.3e-8

// TODO: Perl and awk get a bit ugly in nextflow, move to a proper external script
process execute {
    tag "${name}"
    input:
        tuple val(name), path(harvest), path(finder), path(genome)
    output:
        tuple val(name), path("${genome.baseName}.pass.list"), optional: true
    conda 'bioconda::ltr_retriever'
    cpus 8
    memory 12.GB
    time '8h'
    publishDir 'out_ltr_retriever'
shell:
"""
cat ${harvest} ${finder} > rawLTR.scn
LTR_retriever -name ${name} -inharvest rawLTR.scn \
    -genome ${genome} \
    -u ${params.u} \
    -threads ${task.cpus} \
    -noanno
"""
}

process process_output {
    tag "${name}"
    input:
        tuple val(name), path(harvest)
    output:
        tuple val(name), path("${name}.LTR.intact.raw.fa", optional: true), path("${name}.LTRlib.fa", optional: true), path("${name}.LTR.intact.raw.fa.anno.list", optional: true)
    conda 'bioconda::trf bioconda::mdust bioconda::tesorter'
    cpus 1
    memory 4.GB
    time '1h'
    publishDir 'out_ltr_retriever'
shell:
'''
awk '{if (\$1 !~ /#/) print \$1\"\\t\"\$1}' !{name}.pass.list |\
perl !{projectDir}/util/call_seq_by_list.pl - -C !{name} > !{name.baseName}.LTR.intact.fa.ori
perl -i -nle 's/\\|.*//; print \$_' !{name.baseName}.LTR.intact.fa.ori
perl !{projectDir}/util/rename_LTR_skim.pl !{name.baseName}.LTR.intact.fa.ori !{name}.defalse > !{name.baseName}.LTR.intact.fa.anno
mv !{name.baseName}.LTR.intact.fa.anno !{name.baseName}.LTR.intact.fa.ori

# remove simple repeats and candidates with simple repeats at terminals
mdust !{name.baseName}.LTR.intact.fa.ori > !{name.baseName}.LTR.intact.fa.ori.dusted
perl !{projectDir}/util/cleanup_tandem.pl -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 100 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f !{name.baseName}.LTR.intact.fa.ori.dusted > !{name.baseName}.LTR.intact.fa.ori.dusted.cln

# annotate and remove not LTR candidates
TEsorter !{name.baseName}.LTR.intact.fa.ori.dusted.cln --disable-pass2 -p !{task.cpus}
perl !{projectDir}/util/cleanup_misclas.pl !{name.baseName}.LTR.intact.fa.ori.dusted.cln.rexdb.cls.tsv
mv !{name.baseName}.LTR.intact.fa.ori.dusted.cln.cln !{name.baseName}.LTR.intact.raw.fa
mv !{name.baseName}.LTR.intact.fa.ori.dusted.cln.cln.list !{name.baseName}.LTR.intact.raw.fa.anno.list

# generate annotated output and gff
perl !{projectDir}/util/output_by_list.pl 1 \
    !{name.baseName}.LTR.intact.fa.ori 1 \
    !{name.baseName}.LTR.intact.raw.fa \
    -FA -ex|grep \\>|perl -nle \
    's/>//; print "Name\\t\$_"' > !{name.baseName}.LTR.intact.fa.ori.rmlist

perl !{projectDir}/util/filter_gff3.pl !{name.baseName}.pass.list.gff3 \
    !{name.baseName}.LTR.intact.fa.ori.rmlist \
    | perl -nle 's/LTR_retriever/EDTA/gi; print \$_' \
    > !{name.baseName}.LTR.intact.raw.gff3
'''
}

workflow LTR_RETRIEVER {
    take:
        input
    
    main:
        input | execute | process_output

    emit:
        process_output.out
        
}
