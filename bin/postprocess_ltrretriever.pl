#!/usr/bin/env perl

use warnings;
use strict;

my $genome_arg = $ARGV[0];
my $threads = $ARGV[1];
my $bin_path = $ARGV[2];

# Assumptions
my $genome = $genome_arg;
my $call_seq = "$bin_path/call_seq_by_list.pl";
my $rename_LTR = "$bin_path/rename_LTR_skim.pl";
my $mdust = '';
my $cleanup_tandem = "$bin_path/cleanup_tandem.pl";
my $trf = "`which trf`"; # $trf = '' not working with arg parsing logic
my $TEsorter = '';
my $cleanup_misclas = "$bin_path/cleanup_misclas.pl";
my $output_by_list = "$bin_path/output_by_list.pl";
my $filter_gff = "$bin_path/filter_gff3.pl";

# Take
# Start:    https://github.com/jguhlin/EDTA/blob/59d28bae71061c1a4fd5c5c848d9915bb131dd99/EDTA_raw.pl#L366
# End:      https://github.com/jguhlin/EDTA/blob/59d28bae71061c1a4fd5c5c848d9915bb131dd99/EDTA_raw.pl#L389

# Changes
# - Removed multiline backticks
# - Commented out cp $genome.LTR.intact.raw.fa.anno.list ../
# - COmmented out `rm $genome`


# get full-length LTR from pass.list
`awk '{if (\$1 !~ /#/) print \$1"\\t"\$1}' $genome.pass.list | perl $call_seq - -C $genome > $genome.LTR.intact.fa.ori`;
`perl -i -nle 's/\\|.*//; print \$_' $genome.LTR.intact.fa.ori`;
`perl $rename_LTR $genome.LTR.intact.fa.ori $genome.defalse > $genome.LTR.intact.fa.anno`;
`mv $genome.LTR.intact.fa.anno $genome.LTR.intact.fa.ori`;

# remove simple repeats and candidates with simple repeats at terminals
`${mdust}mdust $genome.LTR.intact.fa.ori > $genome.LTR.intact.fa.ori.dusted`;
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.9 -minlen 100 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.LTR.intact.fa.ori.dusted > $genome.LTR.intact.fa.ori.dusted.cln`;

# annotate and remove not LTR candidates
`${TEsorter}TEsorter $genome.LTR.intact.fa.ori.dusted.cln --disable-pass2 -p $threads 2>/dev/null`;
`perl $cleanup_misclas $genome.LTR.intact.fa.ori.dusted.cln.rexdb.cls.tsv`;
`mv $genome.LTR.intact.fa.ori.dusted.cln.cln $genome.LTR.intact.raw.fa`;
`mv $genome.LTR.intact.fa.ori.dusted.cln.cln.list $genome.LTR.intact.raw.fa.anno.list`;
# `cp $genome.LTR.intact.raw.fa.anno.list ../`;

# # generate annotated output and gff
`perl $output_by_list 1 $genome.LTR.intact.fa.ori 1 $genome.LTR.intact.raw.fa -FA -ex|grep \\>|perl -nle 's/>//; print "Name\\t\$_"' > $genome.LTR.intact.fa.ori.rmlist`;
`perl $filter_gff $genome.pass.list.gff3 $genome.LTR.intact.fa.ori.rmlist | perl -nle 's/LTR_retriever/EDTA/gi; print \$_' > $genome.LTR.intact.raw.gff3`;
# `rm $genome`
