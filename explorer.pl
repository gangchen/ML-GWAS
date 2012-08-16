
use warnings;
use strict;


open(my $file, "psoriasis_all.snp") or die $!;
#open(my $file, "test.data") or die $!;
my %strings;
my $linenum = 0;
while(my $line = <$file>){
    chomp($line);
    my ($chr, $pos, $ref, $label, @samples) = split "\t",$line;
    my @characters;
    my %num;
    
    for(@samples){
        if(m/(.)(\d+)(.)(\d+)(.)(\d+)/g){
            my $base = $1;
            
            my $quality = $2;
            my $deepth = $4+$6;
            if($quality >= 20 and $deepth >= 4){
                push @characters, $1;    
                $num{$base}++;
            }else{
                push @characters, "0";
            }
	   }
    }
    
    my %map;
    my $num = 1;
    $map{"0"} = 0;
    for(keys %num){
	   $map{$_} = $num++;
    }
    next if($num <= 2); # check if this site is a SNP
    print $chr."-".$pos;
    for(@characters){
	print "\t", $map{$_};
    }
    print "\n";
}
close $file;
