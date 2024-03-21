// Output files
// !{name.baseName}.LTR.intact.raw.fa
// !{name.baseName}.LTRlib.fa
// !{name.baseName}.LTR.intact.fa // seems to be optional
// !{name.baseName}.LTR.intact.gff3 // seems to be optional
// !{name.baseName}.LTR.intact.raw.fa.anno.list

params.u                = 1.3e-8

// TODO: Perl and awk get a bit ugly in nextflow, move to a proper external script
process ltr_retriever {
    tag "${name}"
    input:
        tuple val(name), path(harvest), path(finder), path(genome)
    output:
        tuple val(name), path("${genome}.pass.list"), path("${genome}.defalse"), path("${genome}.pass.list.gff3"), path("${genome}.LTRlib.fa"), path(genome), optional: true
    conda 'bioconda::ltr_retriever'
    cpus 4
    memory 8.GB
    time '4h'
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
        tuple val(name), path(retriever), path(defalse), path(passlist), path(ltrlib), path(assembly)
    output:
        tuple val(name), path("${assembly}.LTR.raw.fa"), path("${assembly}.LTR.intact.raw.gff3"), path("${assembly}.LTR.intact.raw.fa"), path("${assembly}.LTR.intact.gff3"), path("${assembly}.LTR.intact.gff3"), optional: true
    conda 'bioconda::trf bioconda::mdust bioconda::tesorter'
    cpus 1
    memory 4.GB
    time '1h'
    publishDir 'out_ltr_retriever'
shell:
'''
awk '{if (\$1 !~ /#/) print \$1\"\\t\"\$1}' !{retriever} |\
perl !{projectDir}/util/call_seq_by_list.pl - -C !{assembly} > !{name}.LTR.intact.fa.ori
perl -i -nle 's/\\|.*//; print \$_' !{name}.LTR.intact.fa.ori
perl !{projectDir}/util/rename_LTR_skim.pl !{name}.LTR.intact.fa.ori !{defalse} > !{name}.LTR.intact.fa.anno
mv !{name}.LTR.intact.fa.anno !{name}.LTR.intact.fa.ori

# remove simple repeats and candidates with simple repeats at terminals
mdust !{name}.LTR.intact.fa.ori > !{name}.LTR.intact.fa.ori.dusted
perl !{projectDir}/util/cleanup_tandem.pl -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 100 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f !{name}.LTR.intact.fa.ori.dusted > !{name}.LTR.intact.fa.ori.dusted.cln

# annotate and remove not LTR candidates
TEsorter !{name}.LTR.intact.fa.ori.dusted.cln --disable-pass2 -p !{task.cpus}
perl !{projectDir}/util/cleanup_misclas.pl !{name}.LTR.intact.fa.ori.dusted.cln.rexdb.cls.tsv
mv !{name}.LTR.intact.fa.ori.dusted.cln.cln !{name}.LTR.intact.raw.fa
mv !{name}.LTR.intact.fa.ori.dusted.cln.cln.list !{name}.LTR.intact.raw.fa.anno.list

# generate annotated output and gff
perl !{projectDir}/util/output_by_list.pl 1 \
    !{name}.LTR.intact.fa.ori 1 \
    !{name}.LTR.intact.raw.fa \
    -FA -ex|grep \\>|perl -nle \
    's/>//; print "Name\\t\$_"' > !{name}.LTR.intact.fa.ori.rmlist

perl !{projectDir}/util/filter_gff3.pl !{passlist} \
    !{name}.LTR.intact.fa.ori.rmlist \
    | perl -nle 's/LTR_retriever/EDTA/gi; print \$_' \
    > !{name}.LTR.intact.raw.gff3

cp !{ltrlib} !{name}.LTR.raw.fa
'''
}

workflow LTR_RETRIEVER {
    take:
        input
    
    main:
        input | ltr_retriever | process_output

    emit:
        process_output.out
        
}
