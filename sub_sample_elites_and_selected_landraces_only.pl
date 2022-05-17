#!/usr/bin/perl

open(IN, "core_watkins.txt");
while(<IN>){chomp; $core{$_}++;}
$n = keys %core;
close IN;


open(IN, "high_heterozygosity_lines.txt");
while(<IN>){chomp; ($line, $hetrate) = split(/\t/,$_); $avoid{$line}++;}
$n = keys %avoid;
close IN;


$chr = "$ARGV[0]";
chomp $chr;

#Substitute in the path to your own VCF files here.  
open(IN, "/home/bzglab/rdsf/elite_plus_watseq_raw_vcf/$chr.vcf" ||die "Can't find any vcf file data in the path specified.");
`mkdir $chr`;

open(VCFHEAD, ">$chr/$chr.vcf.header"); print VCFHEAD "##fileformat=VCFv4.2\n";

open(OUT, ">$chr/${chr}.vcf");

open(LOG, ">$chr/$chr.sub-sample.log");
while(<IN>)
{


if(/^#CHROM/)
 {
 $start = 1;
 chomp;
 ($id, @header) = split(/\t/, $_);
 print OUT "$id"; print VCFHEAD "$id";
 $i = 0;
 foreach $cell (@header)
   {
   if( ($cell !~ /^WATDE\d/ || $core{$cell}>0) && $avoid{$cell} <1){print OUT "\t$cell"; print VCFHEAD "\t$cell"; $non_watkins++; $selected_cols{$i}++;}
   $i++;
   }
 print OUT "\n";
 print VCFHEAD "\n";
 $n = keys %selected_cols;
 #print "There are $n selected columns\n"; 
 }

if($start <1 ){print OUT;}

if($start >0 && $_ !~ /^#CHROM/)
 {
 $l++;
 $i = 0;
 chomp;
 ($id, @data) = split(/\t/, $_);
 print OUT "$id";
 foreach $cell (@data){if($selected_cols{$i}>0){print OUT "\t$cell"} $i++;}
 print OUT "\n";
 }

}
print LOG "Data written for $non_watkins varieties and $l SNP lines\n";

close OUT;


