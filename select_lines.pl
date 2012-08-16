use warnings;
use strict;

my @file;
open(my $file, "4X_20_80_2bit_psoriasis.txt") or die $!;
my @content = <$file>;
close $file;

open(my $anno, "4X_20_80_psoriasis_annovar.txt.exonic_variant_function") or die $!;
open(my $output, ">4X_20_80_2bit_psoriasis_nonsynonymous.txt") or die $!;
while(my $line = <$anno>){
	my ($line, $label, $genename, $chr, $pos, @other) = split "\t", $line;
	if($label eq "nonsynonymous SNV"){
		$line =~ s/line//;
		$genename =~ s/(.+?):.+/$1/;
		print $output $genename, "\t",  $content[$line*2 -2 ];
		print $output $genename, "\t",  $content[$line*2 -1 ];
	}
}
close $output;
close $anno;