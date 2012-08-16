use warnings;
use strict;


open(my $refFile, "refGene.txt") or die $!;
my %ref;
while(my $line = <$refFile>){
	my @fields = split "\t", $line;
	my $chr = $fields[2];
	my $start = $fields[4];
	my $end = $fields[5];
	my $gene = $fields[12];
	#print $fields[2], "\t", $fields[4], "\t", $fields[5], "\n";
	$ref{$chr}{$gene}{"start"} = $start;
	$ref{$chr}{$gene}{"end"} = $end;
}

close $refFile;

open(my $snpFile, "4X_20_psoriasis.full.data") or die $!;
while(my $line = <$snpFile>){
	chomp($line);
	my ($snpName, @sample) = split " ",$line;
	my ($chr, $pos) = split "-", $snpName;
	for my $gene (keys %{$ref{$chr}}){
		if($pos >= $ref{$chr}{$gene}{"start"} and $pos <= $ref{$chr}{$gene}{"end"}){
			print $snpName, "\t", $gene, "\t", join("\t", @sample), "\n";
			last;
		}
	}
}

close($snpFile);