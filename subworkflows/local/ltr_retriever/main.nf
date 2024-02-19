// Output files
// !{genome.baseName}.LTR.intact.raw.fa
// !{genome.baseName}.LTRlib.fa
// !{genome.baseName}.LTR.intact.fa // seems to be optional
// !{genome.baseName}.LTR.intact.gff3 // seems to be optional
// !{genome.baseName}.LTR.intact.raw.fa.anno.list

params.u                = 1.3e-8

// TODO: Perl and awk get a bit ugly in nextflow, move to a proper external script
process execute {
    tag "${genome.baseName}"
    input:
        path(genome)
        path(harvest)
        path(finder)
    output:
        tuple path(genome), path ("${genome}.pass.list"), optional: true
    conda 'bioconda::ltr_retriever'
    cpus 8
    memory 12.GB
    time '8h'
    publishDir 'out_ltr_retriever'
shell:
"""
cat ${harvest} ${finder} > rawLTR.scn
LTR_retriever -genome ${genome} -inharvest rawLTR.scn \
    -u ${params.u} \
    -threads ${task.cpus} \
    -noanno
"""
}

process process_output {
    tag "${genome.baseName}"
    input:
        tuple path(genome), path(harvest)
    output:
        tuple path("${genome.baseName}.LTR.intact.raw.fa", optional: true), path("${genome.baseName}.LTRlib.fa", optional: true), path("${genome.baseName}.LTR.intact.raw.fa.anno.list", optional: true)
    conda 'bioconda::trf bioconda::mdust bioconda::tesorter'
    cpus 1
    memory 4.GB
    time '1h'
    publishDir 'out_ltr_retriever'
shell:
'''
awk '{if (\$1 !~ /#/) print \$1\"\\t\"\$1}' !{genome}.pass.list |\
perl !{projectDir}/util/call_seq_by_list.pl - -C !{genome} > !{genome.baseName}.LTR.intact.fa.ori
perl -i -nle 's/\\|.*//; print \$_' !{genome.baseName}.LTR.intact.fa.ori
perl !{projectDir}/util/rename_LTR_skim.pl !{genome.baseName}.LTR.intact.fa.ori !{genome}.defalse > !{genome.baseName}.LTR.intact.fa.anno
mv !{genome.baseName}.LTR.intact.fa.anno !{genome.baseName}.LTR.intact.fa.ori

# remove simple repeats and candidates with simple repeats at terminals
mdust !{genome.baseName}.LTR.intact.fa.ori > !{genome.baseName}.LTR.intact.fa.ori.dusted
perl !{projectDir}/util/cleanup_tandem.pl -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 100 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f !{genome.baseName}.LTR.intact.fa.ori.dusted > !{genome.baseName}.LTR.intact.fa.ori.dusted.cln

# annotate and remove not LTR candidates
TEsorter !{genome.baseName}.LTR.intact.fa.ori.dusted.cln --disable-pass2 -p !{task.cpus}
perl !{projectDir}/util/cleanup_misclas.pl !{genome.baseName}.LTR.intact.fa.ori.dusted.cln.rexdb.cls.tsv
mv !{genome.baseName}.LTR.intact.fa.ori.dusted.cln.cln !{genome.baseName}.LTR.intact.raw.fa
mv !{genome.baseName}.LTR.intact.fa.ori.dusted.cln.cln.list !{genome.baseName}.LTR.intact.raw.fa.anno.list

# generate annotated output and gff
perl !{projectDir}/util/output_by_list.pl 1 \
    !{genome.baseName}.LTR.intact.fa.ori 1 \
    !{genome.baseName}.LTR.intact.raw.fa \
    -FA -ex|grep \\>|perl -nle \
    's/>//; print "Name\\t\$_"' > !{genome.baseName}.LTR.intact.fa.ori.rmlist

perl !{projectDir}/util/filter_gff3.pl !{genome.baseName}.pass.list.gff3 \
    !{genome.baseName}.LTR.intact.fa.ori.rmlist \
    | perl -nle 's/LTR_retriever/EDTA/gi; print \$_' \
    > !{genome.baseName}.LTR.intact.raw.gff3

'''
}

workflow LTR_RETRIEVER {
    take:
        genome
        ltr_harvest_out
        ltr_finder_out
    
    main:
        execute(genome, ltr_harvest_out, ltr_finder_out)
        
        // Filter out
        process_output(ltr_retriever.out, genome)
        
}
