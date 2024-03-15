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
        tuple val("${x.baseName}"), path("${x.baseName}_sanitized.fasta"), path("${x.baseName}_sanitized.translation_table.tsv")
    publishDir 'sanitized_genomes'
    time "10m"
    memory 3.GB
    cpus 1

"""
fffx length-filter ${x} filtered.fa 1000
fffx sanitize filtered.fa ${x.baseName}_sanitized
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
include { ANNOSINE      } from './subworkflows/local/annosine'
include { REPEATMODELER } from './subworkflows/local/repeatmodeler'
include { TIRLEARNER    } from './subworkflows/local/tir_learner'

workflow {
    genomes = channel.fromPath(params.genomes + "/*")
    sanitize(genomes)
 

    // Handy for data control
    genomes_all = sanitize.out.map({ it ->
        tuple(it[0], ["name": it[0], "assembly": it[1], "translation_table": it[2]], it[1])
    })

    genomes_input = sanitize.out.map({ it ->
        tuple(["name": it[0], "assembly": it[1], "translation_table": it[2]], it[1])
    })

    genomes_files = sanitize.out.map({ it ->
        tuple(it[0], it[1])
    })
    
    // Input for (nearly all) functions are going to be the map "data"
    // Output will be a tuple (data.name, output_file)
    // So we can merge them all later using data.name

    // The map has the folliwing features:
    // - name: The name of the genome
    // - assembly: The sanitized genome
    // - translation_table: The translation table for the genome

    // These two feed into ltr_retriever (but the first two run in parallel)
    LTR_HARVEST(genomes_input)
    LTR_FINDER(genomes_input)
    LTR_HARVEST.out.join(LTR_FINDER.out).join(genomes_files).set { ltr_retriever_input }
    LTR_RETRIEVER(ltr_retriever_input)

    // These can also run in parallel
    ANNOSINE(genomes_input)
    REPEATMODELER(genomes_input)
    TIRLEARNER(genomes_input)
    // helitron_scanner(sanitized_genomes)
}
