while (<STDIN>) {
	chomp;
	$_ =~ s/^(>.*)\s+(.*)$/$1#DNA\/Helitron\t\$2/;
	print $_ . "\n";
}
