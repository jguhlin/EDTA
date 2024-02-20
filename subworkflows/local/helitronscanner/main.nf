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

workflow HELITRONSCANNER {
    take:
        path(genome)
    
}