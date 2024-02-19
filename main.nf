#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.genomes          = 'genomes'
params.species          = 'others'
params.cds              = ''
params.curatedlib       = ''
params.rmlib            = ''
params.sensitive        = false
params.anno             = false
params.rmout            = ''
params.maxdiv           = 40
params.evaluate         = true
params.exclude          = ''
params.maxint           = 5000

// TODO: Check parameters, also help
/*

if ($maxdiv < 0 or $maxdiv > 100){die "The expected value for the div parameter is 0 - 100!\n"}
if ($overwrite != 0 and $overwrite != 1){ die "The expected value for the overwrite parameter is 0 or 1!\n"}
if ($sensitive != 0 and $sensitive != 1){ die "The expected value for the sensitive parameter is 0 or 1!\n"}
if ($anno != 0 and $anno != 1){ die "The expected value for the anno parameter is 0 or 1!\n"}
if ($evaluate != 0 and $evaluate != 1){ die "The expected value for the evaluate parameter is 0 or 1!\n"}
if ($force != 0 and $force != 1){ die "The expected value for the force parameter is 0 or 1!\n"}
if ($miu !~ /[0-9\.e\-]+/){ die "The expected value for the u parameter is float value without units!\n"}
if ($debug != 0 and $debug != 1){ die "The expected value for the debug parameter is 0 or 1!\n"}
if ($threads !~ /^[0-9]+$/){ die "The expected value for the threads parameter is an integer!\n"}

*/

// TODO: Check inputed repeat libraries, CDS, etc...
// TODO: Check exclude file

// Rename FASTA headers (just makes everything easier later)
// TODO: Put fffx on bioconda or somewhere so it just runs, otherwise tiny container
// NEW: Filter out small contigs (1000bp or less, remove)
// TODO: Make sanitize command do filtering, bbmap is memory hungry!
process sanitize {
    tag "${x.baseName}"
    input:
        path x
    output:
        tuple path("${x.baseName}_sanitized.fasta"), path("${x.baseName}_sanitized.translation_table.tsv")
    publishDir 'sanitized_genomes'
    time "10m"
    memory 3.GB
    cpus 1

"""
fffx length-filter ${x} filtered.fa 1000
fffx sanitize filtered.fa ${x.baseName}_sanitized
"""
}

process annosine {
    tag "${genome.baseName}"
    input:
        path(genome)
    output:
        path("${genome.baseName}.SINE.raw.fa")
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
  -auto 1 3 !{genome} .

  if [ -s "Seed_SINE.fa" ]; then
    awk '{print \$1}' Seed_SINE.fa > !{genome.baseName}.AnnoSINE.raw.fa
    TEsorter !{genome.baseName}.AnnoSINE.raw.fa --disable-pass2 -p !{task.cpus}
    touch !{genome.baseName}.AnnoSINE.raw.fa.rexdb.cls.tsv
    perl !{projectDir}/util/cleanup_misclas.pl !{genome.baseName}.AnnoSINE.raw.fa.rexdb.cls.tsv

    perl !{projectDir}/util/cleanup_tandem.pl -misschar N \
        -nc 50000 \
        -nr 0.9 \
        -minlen 100 \
        -minscore 3000 \
        -trf 1 \
        -cleanN 1 \
        -cleanT 1 \
        -f !{genome.baseName}.AnnoSINE.raw.fa.cln > !{genome.baseName}.SINE.raw.fa
    else
        touch !{genome.baseName}.SINE.raw.fa
    fi

'''
}

process repeatmodeler {
    tag "${genome.baseName}"
    input:
        path(genome)
    output:
        tuple path("${genome.baseName}.LINE.raw.fa"), path("${genome.baseName}.RM2.fa"), path("${genome.baseName}.RM2.raw.fa")
    conda 'bioconda::repeatmodeler=2.0.2a=pl5262h9ee0642_0 tesorter=1.4.6=pyhdfd78af_0'
    cpus 32
    time '6d'
    memory 32.GB
    publishDir 'out_repeatmodeler'

shell:
'''
BuildDatabase -name !{genome.baseName} !{genome}
RepeatModeler -engine ncbi -threads !{task.cpus} -database !{genome.baseName}

awk '{print \$1}' !{genome.baseName}-families.fa > !{genome.baseName}.RM2.raw.fa
TEsorter !{genome.baseName}.RM2.raw.fa --disable-pass2 -p !{task.cpus}

perl !{projectDir}/util/cleanup_misclas.pl !{genome.baseName}.RM2.raw.fa.rexdb.cls.tsv

TEsorter !{genome.baseName}.RM2.raw.fa.cln --disable-pass2 -p !{task.cpus}
perl -nle 's/>\\S+\\s+/>/; print \$_' !{genome.baseName}.RM2.raw.fa.cln.rexdb.cls.lib > !{genome.baseName}.RM2.raw.fa.cln

perl !{projectDir}/util/cleanup_tandem.pl -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 80 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f !{genome.baseName}.RM2.raw.fa.cln > !{genome.baseName}.RM2.raw.fa.cln2

grep -P 'LINE|SINE' !{genome.baseName}.RM2.raw.fa.cln2 | perl !{projectDir}/util/output_by_list.pl 1 !{genome.baseName}.RM2.raw.fa.cln2 1 - -FA > !{genome.baseName}.LINE.raw.fa
grep -P 'LINE|SINE' !{genome.baseName}.RM2.raw.fa.cln2 | perl !{projectDir}/util/output_by_list.pl 1 !{genome.baseName}.RM2.raw.fa.cln2 1 - -FA -ex > !{genome.baseName}.RM2.fa
'''
}

process tir_learner {
    tag "${genome.baseName}"
    input:
        path(genome)
    output:
        tuple path("${genome.baseName}.TIR.intact.fa"), path("${genome.baseName}.TIR.intact.raw.fa"), path("${genome.baseName}.TIR.intact.raw.gff3")
    cpus 16
    memory 16.GB
    time '12h'
    conda 'python=3.10 bioconda::genometools-genometools bioconda::genericrepeatfinder bioconda::cd-hit bioconda::mdust bioconda::tesorter blast pandas swifter regex scikit-learn tensorflow bioconda::trf'
    publishDir 'out_tir_learner'

"""
python3 ${projectDir}/bin/TIR-Learner3.0/TIR-Learner3.0.py \
    -f ${genome} \
    -s ${params.species} \
    -t ${task.cpus} \
    -l ${params.maxint} \
    -o ${genome.baseName}.TIR

cd ${genome.baseName}.TIR

perl ${projectDir}/util/rename_tirlearner.pl \
    ./TIR-Learner-Result/TIR-Learner_FinalAnn.fa | perl -nle 's/TIR-Learner_//g; print \$_' > ${genome.baseName}.TIR

perl ${projectDir}/util/get_ext_seq.pl ../${genome} ${genome.baseName}.TIR
perl ${projectDir}/util/flanking_filter.pl \
    -genome ../${genome} \
    -query ${genome.baseName}.TIR.ext30.fa \
    -miniden 90 \
    -mincov 0.9 \
    -maxct 20 \
    -t ${task.cpus} \

perl ${projectDir}/util/output_by_list.pl 1 \
    ${genome.baseName}.TIR \
    1 \
    ${genome.baseName}.TIR.ext30.fa.pass.fa \
    -FA -MSU0 -MSU1 \
    > ${genome.baseName}.TIR.ext30.fa.pass.fa.ori

# remove simple repeats and candidates with simple repeats at terminals
mdust ${genome.baseName}.TIR.ext30.fa.pass.fa.ori > ${genome.baseName}.TIR.ext30.fa.pass.fa.dusted
perl ${projectDir}/util/cleanup_tandem.pl \
    -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 80 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f ${genome.baseName}.TIR.ext30.fa.pass.fa.dusted > ${genome.baseName}.TIR.ext30.fa.pass.fa.dusted.cln

# annotate and remove non-TIR candidates
TEsorter ${genome.baseName}.TIR.ext30.fa.pass.fa.dusted.cln --disable-pass2 -p ${task.cpus}
perl ${projectDir}/util/cleanup_misclas.pl \
    ${genome.baseName}.TIR.ext30.fa.pass.fa.dusted.cln.rexdb.cls.tsv
mv ${genome.baseName}.TIR.ext30.fa.pass.fa.dusted.cln.cln ../${genome.baseName}.TIR.intact.raw.fa
cp ${genome.baseName}.TIR.ext30.fa.pass.fa.dusted.cln.cln.list ../${genome.baseName}.TIR.intact.raw.fa.anno.list

# get gff3 of intact TIR elements
perl -nle 's/\\\\-\\\\+\\\\-/_Len:/; my (\$chr, \$method, \$supfam, \$s, \$e, \$anno) = (split)[0,1,2,3,4,8]; my \$class=\"DNA\"; \$class=\"MITE\" if \$e-\$s+1 <= 600; my (\$tir, \$iden, \$tsd)=(\$1, \$2/100, \$3) if \$anno=~/TIR:(.*)_([0-9.]+)_TSD:([a-z0-9._]+)_LEN/i; print \"\$chr \$s \$e \$chr:\$s..\$e \$class/\$supfam structural \$iden . . . TSD=\$tsd;TIR=\$tir\"' ./TIR-Learner-Result/TIR-Learner_FinalAnn.gff3 |\
    perl ${projectDir}/util/output_by_list.pl 4 - 1 \
    ../${genome.baseName}.TIR.intact.raw.fa -MSU0 -MSU1 \
    > ${genome.baseName}.TIR.intact.raw.bed

perl ${projectDir}/util/bed2gff.pl ${genome.baseName}.TIR.intact.raw.bed TIR > ../${genome.baseName}.TIR.intact.raw.gff3
cp ../${genome.baseName}.TIR.intact.raw.fa ../${genome.baseName}.TIR.intact.fa
"""
}

// It's actually single threaded
process helitron_scanner {
    tag "${genome.baseName}"
    input:
        path(genome)
    output:
        path("${genome.baseName}.Helitron.intact.raw.gff3")
    conda 'bioconda::tesorter bioconda::mdust bioconda::trf'
    cpus 1 
    time '18h'
    memory 32.GB
    publishDir 'out_helitron_scanner'

"""
sh ${projectDir}/util/run_helitron_scanner.sh \
    ${genome} \
    ${task.cpus}

perl ${projectDir}/util/format_helitronscanner_out.pl \
    -genome $genome \
    -sitefilter 1 \
    -minscore 12 \
    -keepshorter 1 \
    -extlen 30 \
    -extout 1

perl ${projectDir}/util/flanking_filter.pl \
    -genome ${genome} \
    -query ${genome}.HelitronScanner.filtered.ext.fa \
    -miniden 90 \
    -mincov 0.9 \
    -maxct 5 \
    -t ${task.cpus}

# remove simple repeats and candidates with simple repeats at terminals
perl ${projectDir}/util/output_by_list.pl 1 \
    ${genome}.HelitronScanner.filtered.fa \
    1 \
    ${genome}.HelitronScanner.filtered.ext.fa \
    -FA > ${genome}.HelitronScanner.filtered.pass.fa

mdust ${genome}.HelitronScanner.filtered.pass.fa > ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted
perl ${projectDir}/util/cleanup_tandem.pl \
    -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 100 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted | perl ${projectDir}/util/helitron_renamer.pl > ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln
# Too complicated for nextflow
# perl -nle 's/^(>.*)\\s+(.*)\\\$/\\\$1#DNA\\/Helitron\\t\\\$2/; print \\\$_' > ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln

# annotate and remove non-Helitron candidates
TEsorter ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln \
    --disable-pass2 \
    -p ${task.cpus}

perl ${projectDir}/util/cleanup_misclas.pl \
    ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln.rexdb.cls.tsv
mv ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln.cln ${genome.baseName}.Helitron.intact.raw.fa
cp ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln.cln.list ${genome.baseName}.Helitron.intact.raw.fa.anno.list

# get intact Helitrons and gff3
perl ${projectDir}/util/make_bed_with_intact.pl \
    ${genome.baseName}.Helitron.intact.raw.fa > ${genome.baseName}.Helitron.intact.raw.bed
perl ${projectDir}/util/bed2gff.pl ${genome.baseName}.Helitron.intact.raw.bed HEL > ${genome.baseName}.Helitron.intact.raw.gff3
"""
}

// Includes
include { LTR_HARVEST   } from './subworkflows/local/ltr_harvest'
include { LTR_FINDER    } from './subworkflows/local/ltr_finder'
include { LTR_RETRIEVER } from './subworkflows/local/ltr_retriever'

workflow {
    genomes = channel.fromPath(params.genomes + "/*")
    sanitize(genomes)
 
    // Get first element of each tuple
    sanitized_genomes = sanitize.out.map({ it[0] })

    // These two feed into ltr_retriever (but the first two run in parallel)
    LTR_HARVEST(sanitized_genomes)
    LTR_FINDER(sanitized_genomes)
    LTR_RETRIEVER(sanitized_genomes, LTR_HARVEST.out, LTR_FINDER.out)

    // These can also run in parallel
    annosine(sanitized_genomes)
    repeatmodeler(sanitized_genomes)
    tir_learner(sanitized_genomes)
    helitron_scanner(sanitized_genomes)
}
