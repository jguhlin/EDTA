#!/usr/bin/env perl

use warnings;
use strict;

my $genome_arg = $ARGV[0];
my $threads = $ARGV[1];
my $bin_path = $ARGV[2];

# Assumptions
my $genome = $genome_arg;
my $cleanup_tandem = "$bin_path/cleanup_tandem.pl";
my $trf = "`which trf`"; # $trf = '' not working with arg parsing logic
my $TEsorter = '';
my $cleanup_misclas = "$bin_path/cleanup_misclas.pl";
my $output_by_list = "$bin_path/output_by_list.pl";

# Take
# Start:    https://github.com/jguhlin/EDTA/blob/59d28bae71061c1a4fd5c5c848d9915bb131dd99/EDTA_raw.pl#L515
# End:      https://github.com/jguhlin/EDTA/blob/59d28bae71061c1a4fd5c5c848d9915bb131dd99/EDTA_raw.pl#L527


# annotate and remove misclassified candidates
`awk '{print \$1}' $genome-families.fa > $genome.RM2.raw.fa` if -e "$genome-families.fa";
`${TEsorter}TEsorter $genome.RM2.raw.fa --disable-pass2 -p $threads 2>/dev/null`;
`perl $cleanup_misclas $genome.RM2.raw.fa.rexdb.cls.tsv`;

# reclassify clean candidates
`${TEsorter}TEsorter $genome.RM2.raw.fa.cln --disable-pass2 -p $threads 2>/dev/null`;
`perl -nle 's/>\\S+\\s+/>/; print \$_' $genome.RM2.raw.fa.cln.rexdb.cls.lib > $genome.RM2.raw.fa.cln`;

# clean up tandem repeat
`perl $cleanup_tandem -misschar N -nc 50000 -nr 0.8 -minlen 80 -minscore 3000 -trf 1 -trf_path $trf -cleanN 1 -cleanT 1 -f $genome.RM2.raw.fa.cln > $genome.RM2.raw.fa.cln2`;
`grep -P 'LINE|SINE' $genome.RM2.raw.fa.cln2 | perl $output_by_list 1 $genome.RM2.raw.fa.cln2 1 - -FA > $genome.LINE.raw.fa`;
`grep -P 'LINE|SINE' $genome.RM2.raw.fa.cln2 | perl $output_by_list 1 $genome.RM2.raw.fa.cln2 1 - -FA -ex > $genome.RM2.fa`;