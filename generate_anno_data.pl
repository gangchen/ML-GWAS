#!perl

open(my oriDataFile, "psoriasis_all.snp") or die $!;

close $oriDataFile;


open(my $annoFile, "psoriasis_snp_annovar.txt.exonic_variant_function") or die $!;

close $annoFile;

