process tir_learner {
    tag "${data.name}"
    input:
        tuple val(data), path(assembly)
    output:
        tuple val(data), path(assembly), path("${data.name}.TIR")
    cpus 8
    memory 8.GB
    time '6h'
    conda 'python=3.10 bioconda::genometools-genometools bioconda::genericrepeatfinder bioconda::cd-hit bioconda::mdust bioconda::tesorter blast pandas swifter regex scikit-learn tensorflow'

"""
python3 ${projectDir}/bin/TIR-Learner3.0/TIR-Learner3.0.py \
    -f ${assembly} \
    -s ${params.species} \
    -t ${task.cpus} \
    -l ${params.maxint} \
    -o ${data.name}.TIR
"""
}

process process_output {
    tag "${data.name}"
    input:
        tuple val(data), path(assembly), path(tir_directory)
    output:
        tuple val(data.name), path("${data.name}.TIR.intact.fa"), path("${data.name}.TIR.intact.raw.fa"), path("${data.name}.TIR.intact.raw.gff3")
    cpus 2
    memory 2.GB
    time '1h'
    conda 'python=3.10 bioconda::genometools-genometools bioconda::genericrepeatfinder bioconda::cd-hit bioconda::mdust bioconda::tesorter blast pandas regex bioconda::trf'
    publishDir 'out_tir_learner'

"""
# cd ${tir_directory}

perl ${projectDir}/util/rename_tirlearner.pl \
    ${tir_directory}/TIR-Learner-Result/TIR-Learner_FinalAnn.fa |\
    perl -nle 's/TIR-Learner_//g; print \$_' > out.TIR

perl ${projectDir}/util/get_ext_seq.pl ${assembly} out.TIR
perl ${projectDir}/util/flanking_filter.pl \
    -genome ${assembly} \
    -query out.TIR.ext30.fa \
    -miniden 90 \
    -mincov 0.9 \
    -maxct 20 \
    -t ${task.cpus} \

perl ${projectDir}/util/output_by_list.pl 1 \
    out.TIR \
    1 \
    out.TIR.ext30.fa.pass.fa \
    -FA -MSU0 -MSU1 \
    > ${data.name}.TIR.ext30.fa.pass.fa.ori

# remove simple repeats and candidates with simple repeats at terminals
mdust ${data.name}.TIR.ext30.fa.pass.fa.ori > ${data.name}.TIR.ext30.fa.pass.fa.dusted
perl ${projectDir}/util/cleanup_tandem.pl \
    -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 80 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f ${data.name}.TIR.ext30.fa.pass.fa.dusted > ${data.name}.TIR.ext30.fa.pass.fa.dusted.cln

# annotate and remove non-TIR candidates
TEsorter ${data.name}.TIR.ext30.fa.pass.fa.dusted.cln --disable-pass2 -p ${task.cpus}
perl ${projectDir}/util/cleanup_misclas.pl \
    ${data.name}.TIR.ext30.fa.pass.fa.dusted.cln.rexdb.cls.tsv
mv ${data.name}.TIR.ext30.fa.pass.fa.dusted.cln.cln ${data.name}.TIR.intact.raw.fa
cp ${data.name}.TIR.ext30.fa.pass.fa.dusted.cln.cln.list ${data.name}.TIR.intact.raw.fa.anno.list

# get gff3 of intact TIR elements
perl -nle 's/\\\\-\\\\+\\\\-/_Len:/; my (\$chr, \$method, \$supfam, \$s, \$e, \$anno) = (split)[0,1,2,3,4,8]; my \$class=\"DNA\"; \$class=\"MITE\" if \$e-\$s+1 <= 600; my (\$tir, \$iden, \$tsd)=(\$1, \$2/100, \$3) if \$anno=~/TIR:(.*)_([0-9.]+)_TSD:([a-z0-9._]+)_LEN/i; print \"\$chr \$s \$e \$chr:\$s..\$e \$class/\$supfam structural \$iden . . . TSD=\$tsd;TIR=\$tir\"' ${tir_directory}/TIR-Learner-Result/TIR-Learner_FinalAnn.gff3 |\
    perl ${projectDir}/util/output_by_list.pl 4 - 1 \
    ${data.name}.TIR.intact.raw.fa -MSU0 -MSU1 \
    > ${data.name}.TIR.intact.raw.bed

perl ${projectDir}/util/bed2gff.pl ${data.name}.TIR.intact.raw.bed TIR > ${data.name}.TIR.intact.raw.gff3
cp ${data.name}.TIR.intact.raw.fa ${data.name}.TIR.intact.fa
"""
}

workflow TIRLEARNER {
    take:
        genomes

    main:
        genomes | tir_learner | process_output

    emit:
        process_output.out
}