#!/usr/bin/perl

$vcfpath = $ARGV[0];

open(OUT, ">high_heterozygosity_lines.txt");

#Get the line with the variety names - it is the 75th one in our vcf files but might be different in other analyses so check!
$head = `head -75 $vcfpath/chr1A.vcf |tail -1`;
chomp $head;
print "This should be a list of variety names from the CHR1A vcf file:\n$head\n";
if($head !~ /\d/){die "Can't find any vcf file data in the path specified\n";}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	

($chromosome, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @names) = split(/\t/, $head);

#print "First var in column is: $names[0]\n";

@files = `ls */binned_output/out_bin* |grep out_bin_`;

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
