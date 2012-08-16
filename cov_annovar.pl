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
                "W" => "AT","S" => "CG",
                "K" => "TG","M" => "AC",
                "Y" => "TC","R" => "AG",};

# my $encoding = {"aa" => ["1", "0", "0"],
# 				"ab" => ["0", "1", "0"],
# 				"bb" => ["0", "0", "1"],
# 				"0"  => ["0", "0", "0"],};
my $encoding = {"aa" => ["1", "0"],
				"ab" => ["0.5", "0.5"],
				"bb" => ["0", "1"],
				"0"  => ["0", "0"],};


open(my $snpFile, "psoriasis_all.snp") or die $!;
open(my $dataFile, ">4X_20_80_2bit_psoriasis.txt") or die $!;

while(my $line = <$snpFile>){
	chomp $line;
	my @fields = split "\t",$line;
	$fields[0] =~ /chr(.+)/;
	my $chr = $1;
	my $pos = $fields[1];
	my $numCall = 0;
	my %genotypes;
	for my $field ( @fields[4 .. $#fields] ){
		$field =~ /(.)(\d+).(\d+).(\d+)/;
		my $deepth = $3 + $4;
		my $genotype = $1;
		my $quality = $2;
		#print $deepth, "\t", $genotype, "\t", $quality, "\n";
		if($deepth >= 4 and $quality >= 20){
			my ($genotype1, $genotype2) = split "", $map->{$genotype};
			$genotypes{$genotype1}++;
			$genotypes{$genotype2}++;
			$numCall++;
		}
	}

	if(scalar keys %genotypes >= 2 and ($numCall/1500) >= 0.8 ){
		my $majorSNP = "XX";
		my $minorSNP = "XX";
		
		for my $snp (keys %genotypes){
			#print $snp, ":", $genotypes{$snp}, "\t";
			if(!defined($genotypes{$majorSNP}) or 
				$genotypes{$majorSNP} < $genotypes{$snp}){
				$majorSNP = $snp;
			}
		}

		for my $snp (keys %genotypes){
			next if($snp eq $majorSNP);
			if(!defined($genotypes{$minorSNP}) or 
				$genotypes{$minorSNP} < $genotypes{$snp}){
				$minorSNP = $snp;
			}
		}
		 # print $chr, "\t", $pos, "\t", $pos, "\t", $majorSNP, "\t", $minorSNP, "\tcomments: ";
		 # if($genotypes{$minorSNP}/($genotypes{$minorSNP}+$genotypes{$majorSNP}) < 0.01 ){
			# print "rare SNP\n";
		 # }else{
		 # 	print "common SNP\n";
		 # }

		my @features;
		my $feature;
		for my $field ( @fields[4 .. $#fields] ){
			$field =~ /(.)(\d+).(\d+).(\d+)/;
			my $deepth = $3 + $4;
			my $genotype = $1;
			my $quality = $2;
			if($deepth >= 4 and $quality >= 20){
				my ($genotype1, $genotype2) = split "", $map->{$genotype};
				if($genotype1 eq $genotype2){
					if($genotype1 eq $majorSNP){
						push @features, $_ for(@{$encoding->{"aa"}});
					}elsif($genotype1 eq $minorSNP){
						push @features, $_ for(@{$encoding->{"bb"}});
					}else{
						push @features, $_ for(@{$encoding->{"0"}});
					}
				}elsif($genotype1 ne $genotype2 and
						(($genotype1 eq $majorSNP and $genotype2 eq $minorSNP) or
						($genotype2 eq $majorSNP and $genotype1 eq $minorSNP))){
					push @features, $_ for(@{$encoding->{"ab"}});
				}else{
					push @features, $_ for(@{$encoding->{"0"}});
				}
			}else{
				push @features, $_ for(@{$encoding->{"0"}});
			}
		}
		my @nums1 = map {2 * ($_-1)} 1..((scalar @features)/2);
		my @nums2 = map {2 * ($_-1) + 1} 1..((scalar @features)/2);
		#my @nums3 = map {3 * ($_-1) + 2} 1..((scalar @features)/3);
		print $dataFile $chr, "-", $pos, "-1\t", join("\t", @features[@nums1]), "\n";
		print $dataFile $chr, "-", $pos, "-2\t", join("\t", @features[@nums2]), "\n";
		#print $dataFile $chr, "-", $pos, "-3\t", join("\t", @features[@nums3]), "\n";
	}


}
close $dataFile;
close $snpFile;

