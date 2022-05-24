#!/usr/bin/perl


#First load poistions of SNPs from other sources and make sure these are selected at this first stage if they overlap.

($chrfile, @selectedSNPfiles) = @ARGV;


#The format of these csv files should be as per the example in the first three lines shown here:

#Chromosome (IWGSC v1),Position (IWGSC v1),35K_SNPId,Sequence
#chr1A,1174865,AX-94493709,AGGAGTGATGTATCTGAAAATCTCGAGGGAGCTGA[A/G]ACCCAAGCTTCACAAGAACCTCCTGAGTTTCTGAA
#chr1A,1236448,AX-94772289,GCTCGGTCGTCTTCGGCCTTGCGACCAAGCATGAT[C/T]TGGTGGATGGCAGGATGAAGAGGTGGTACCTCAGG



foreach $selectedSNPfile(@selectedSNPfiles)
{
open(IN, "$selectedsnpfile");
while(<IN>)
{
($chr, $pos, @other) = split(/[\t\,]/, $_);
if($chr =~ /^[1-7]/){$chr = "chr$chr";}
$combo = "${chr}_$pos";
#print "$combo\n";
$selected{$combo}++;
}
close IN;
}

$max_het_prop = 0.005;
$min_maf = 0.01;
$min_call_rate = 0.95;

chomp $chrfile;
print "Incoming file name is $file\n";
if($file =~ /(chr[1-7][ABDT]).vcf/){$name = $1; $file = "$name/$name.vcf"; } else{die "Unexpected file name $file\n";}

print "Filtering $file\n";

open(OUT, ">$file.filtered");
open(SUMMARY, ">$name/$name.summary.txt");
print SUMMARY "chr\tid\tqual\thet\thom\tcall1\tcall2\tminor\n";
chomp $file;
open(VCF, "$file" || die "Cant open $file\n");

$printer = 1;
while($line = <VCF>)

{
$i++; $l++;
if($i == 100000){$j++; $i = 0; print "${j}00000 lines processed, $print selected\n";}
chomp $line;
($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format,@data) = split(/\t/, $line);
$combo = "${chr}_$pos";
if($printer ==1){print OUT "$line\n";}
if($line =~ /^\#CHROM/){$printer = 0; $n = @data; $max_het_calls = $n * $max_het_prop; print "processing data for $n varieties\n"; }
if($chr =~ /chr[1-7][ABDT]/)
{
$qual = (int($qual/1000))*1000;
#print "$qual\n";
$quals{$qual}++;
#Qual was set to 50000
if($qual >5000)
{
$het = 0; $hom = 0; $call1 = 0; $call2 = 0; $total =0;
 
foreach $cell(@data)
{
@fields = split(/:/, $cell);
if($fields[0] =~ /^0[\/\|]1/){$het++; if($het >$max_het_calls){last;}}
elsif($fields[0] =~ /^0[\/\|]0/){$hom++; $call1++;}
elsif($fields[0] =~ /^1[\/\|]1/){$hom++; $call2++; }
}

$total = $het+$hom;
if($total > 0)
{
$call_rate = $total/$n;

if($call1 > $call2){$minor = $call2;} else{$minor = $call1;}
if($selected{$combo}<1){$selected{$combo} = 0;}
if($selected{$combo} >0 || ($call_rate >$min_call_rate &&  $het < $max_het_calls && ($minor/$total) >$min_maf))
 {
 print OUT "$line\n";
 print SUMMARY "$chr\t$id\t$qual\t$het\t$hom\t$call1\t$call2\t$minor\t$selected{$combo}\n";
 $print++;
 }
}
}


}
}
print SUMMARY "\n$l SNPs processed\t$print passed filters:\nMax het proportion $max_het_prop\nMin maf $min_maf\nMin call rate $min_call_rate\n";


close OUT; 
close IN;

open(QUAL, ">$name/$name.qual_score_histo.txt");
foreach $qual(sort {$b<=>$a} keys %quals){print QUAL "$qual\t$quals{$qual}\n";}



#Comment out die to proceed and run snpEFF
#$do = 1;

if($do  ==1)
{
#Do snpEff here
$name = $file;
$name =~ s/\.vcf//;

`java -Xmx8g -jar ~/snpEff/snpEff.jar -v -stats $name.html Triticum_aestivum $file.filtered > $name.ann.vcf`;


$print = 0;
open(OUT, ">$name.non_synonymous.vcf"); 
open(VCF, "$name.ann.vcf" || die "Cant open $file\n"); 
$printer = 1; while($line = <VCF>) 
{ 
$i++; if($i == 1000){$j++; $i = 0; print "${j}000 lines processed, $print selected\n";} 
chomp $line; 
($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @data) = split(/\t/, $line); 
if($printer ==1){print OUT "$line\n";} 
if($line =~ /^\#CHROM/){$printer = 0;} if($chr =~ /chr[1-7][ABDT]/ && $info =~ /missense_variant/)
  {
  print OUT "$line\n";
  $print++;
  }
}

print SUMMARY "$print SNPs are misssense\n";
close SUMMARY;
}




#`rm $file.filtered`;
#`rm $file`;

