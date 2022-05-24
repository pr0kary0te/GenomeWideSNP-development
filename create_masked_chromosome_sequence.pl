#!/usr/bin/perl

#This script takes the IWGSC chromosomes and makes a version where each SNP identified will be masked out with N's
# starting $upstream bases upstream and extending $length bases in total.
# The logic is that candidate SNPs and their flanking sequence can then be BWA mapped to this file, and any that have a match must have more than one genome
#match - given that their canonical position has been masked out. This will act as a high copy filter for potential SNPs.


#The chromosome to work on is specified in the command line as a single argument.
$chr = "$ARGV[0]";

#Edit the path to the location of the IWGSC version 1.0 chromosomes.  This directory must be writable   
$path = $ARGV[1];


$upstream = 10;
$length = 21;
$string = "N" x 21;


chomp $chr;


#Get every SNP coordinate for the current chromosome
open(GEN, "$chr/$chr.vcf.filtered.genotypes");
$head = <GEN>;

while(<GEN>)
{
if(/^chr[1-7][ABD]_(\d+)\t/){$positions{$1}++;}
}
close GEN;

$chr =~ s/chr//;

open(CHR, "$path/$chr.fa");
open(OUT, ">$path/$chr.masked.fa");

$head = <CHR>; print OUT "$head";

while(<CHR>)
{
chomp;
$sequence.=$_;
}

close CHR;

foreach $pos(sort {$a<=>$b} keys %positions)
{
$original = substr($sequence, ($pos - $upstream -1), $length);
substr($sequence, ($pos - $upstream -1), $length, $string);
$post = substr($sequence, ($pos - $upstream -1), $length);
#print "$pos $original\n$pos $post\n\n";
} 

print OUT "$sequence\n";

$file = "$path/$chr.masked.fa";
if(-e "$file.sa"){print "$file.sa exists skipping indexing\n"; }
else{ `bwa index $path/$chr.masked.fa`;}

