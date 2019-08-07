#!/usr/bin/perl 

=head
Usage:  perl dose_to_genotype.pl input_csv_file [offset]
        
		Converts MMAP csv file dosages to genotypes (integers)
		offset parameter defaults to 8 (for IMPUTED files with INFO column)
=cut

use strict;

my $datasetLD = $ARGV[0];
my $offset = $ARGV[1];
if($offset == 0) {$offset = 8;}

open(IMPUTED,"<","$datasetLD");
open(GENOTYPES,">","$datasetLD.genotypes");
my $line = <IMPUTED>;
print GENOTYPES $line;  # write header line (no changes)
while($line=<IMPUTED>) {
	chomp $line;
	my @row=split(",",$line);
	my $line_values = @row;
	# print "line-values:[$line_values]\n";
	for (my $i=8; $i < $line_values; $i++) {
		if($row[$i] >= 1.5) {
			$row[$i] = 2;
		} elsif($row[$i] > 0.5) {
			$row[$i] = 1;
		} else { $row[$i] = 0; }
	}
	my $out = join(",",@row);
	print GENOTYPES "$out\n";
}
close(IMPUTED);
close(GENOTYPES);
exit(0);
