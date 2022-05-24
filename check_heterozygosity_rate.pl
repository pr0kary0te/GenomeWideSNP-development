#!/usr/bin/perl

$string = $ARGV[0];
chomp $string;

$vcfpath = $ARGV[1];

open(OUT, ">high_heterozygosity_lines.txt");


if($string =~ /(chr[1-7][ABD])/){$chr = $1;} else{$chr = "*";}


#Substitute in the path to the raw vcf files here
$head = `head -75 $vcfpath/chr1A.vcf |tail -1`;
chomp $head;
print "$head\n";
if($head !~ /\d/){die "Can't find any vcf file data in the path specified\n";}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	

($chromosome, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @names) = split(/\t/, $head);

#print "First var in column is: $names[0]\n";

@files = `ls $chr/binned_output/out_bin* |grep out_bin_`;

foreach $file(@files)
{
chomp $file;
#print "$file\n";
@lines = `cat $file`;
foreach $line(@lines)
 {
 chomp $line;
 ($score, $resolved, $name, $string) = split(/\t/, $line);
 @calls = split(//, $string);
 $i = 0; foreach $call(@calls){$lookup{$i}{$call}++; $i++; $calltypes{$call}++;}
 }
}

print "Var";
foreach $calltype(sort keys %calltypes){print "\t$calltype";}
print "\n";

foreach $i (sort {$a<=>$b} keys %lookup)
{
$var = $names[$i];
print "$var";
foreach $calltype(sort keys %calltypes){print "\t$lookup{$i}{$calltype}";}
$het_rate = $lookup{$i}{1}/($lookup{$i}{0}+$lookup{$i}{2});
print "\t$het_rate\n";
if($het_rate > $max_het_rate){print OUT "$var\t$het_rate\n";}
}
