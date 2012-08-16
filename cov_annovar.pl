=head1 FILE NAME
	cov_annovar.pl

=head1 DESCRIPTION
	script for processing and encoding genotyping data file

=head1 USAGE
	perl cov_annovar.pl
	

=head1 AUTHOR
	Gang Chen, <chengangcs@gmail.com> or <chengang@genomics.cn>
	
=head1 DATE
	2012-08-16
=cut

use warnings;
use strict;			

my $map = {
		"A" => "AA",
		"C" => "CC",
		"G" => "GG",
		"T" => "TT",
        "W" => "AT",
        "S" => "CG",
        "K" => "TG",
        "M" => "AC",
        "Y" => "TC",
        "R" => "AG",};


my $snpFilename = shift;
my $dataFilename = shift;
my $annovarFilename = shift;

open(my $snpFile, $snpFilename) or die $!;
open(my $dataFile, ">".$snpFilename) or die $!;
open(my $annovarFile, ">".$annovarFilename) or die $!;

while(my $line = <$snpFile>){
	chomp $line;
	my @fields = split "\t",$line;
	$fields[0] =~ /chr(.+)/;
	my $chr = $1;
	my $pos = $fields[1];
	my $numCall = 0;
	my %genotypes;
	my @samples;
	for my $field ( @fields[4 .. $#fields] ){
		$field =~ /(.)(\d+)(.)(\d+)(.)(\d+)/;
		my $g1 = $3;
		my $g2 = $5;
		my $d1 = $4;
		my $d2 = $6;
		my %tmp;
		# determine the major allel by number of reads
		my ($genotype1, $genotype2) = split "", $map->{$g1};
		$genotypes{$genotype1} += $d1;
		$genotypes{$genotype2} += $d1;
		$tmp{$genotype1} += $d1;
		$tmp{$genotype2} += $d1;
		($genotype1, $genotype2) = split "", $map->{$g2};
		$genotypes{$genotype1} += $d1;
		$genotypes{$genotype2} += $d1;
		$tmp{$genotype1} += $d2;
		$tmp{$genotype2} += $d2;
		push @samples, \%tmp;
	}

	my $majorSNP = "X";
	my $minorSNP = "X";
	
	for my $snp (keys %genotypes){
		if(!exists($genotypes{$majorSNP}) or 
			$genotypes{$majorSNP} < $genotypes{$snp}){
			$majorSNP = $snp;
		}
	}

	for my $snp (keys %genotypes){
		next if($snp eq $majorSNP);
		if(!exists($genotypes{$minorSNP}) or 
			$genotypes{$minorSNP} < $genotypes{$snp}){
			$minorSNP = $snp;
		}
	}
	print $annovarFile $chr, "\t", $pos, "\t", $pos, "\t", $majorSNP, "\t", $minorSNP, "\tcomments\n";

	print $dataFile $chr, "-", $pos, "-1";
	for(@samples){
		if(exists($_->{$majorSNP})){
			print $dataFile "\t",$_->{$majorSNP};
		}else{
			print $dataFile "\t","0";
		}
	}
	print $dataFile "\n";
	print $dataFile $chr, "-", $pos, "-2";
	for(@samples){
		if(exists($_->{$minorSNP})){
			print $dataFile "\t",$_->{$minorSNP};
		}else{
			print $dataFile "\t","0";
		}
	}
	print $dataFile "\n";

}
close $annovarFile;
close $dataFile;
close $snpFile;

