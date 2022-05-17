#!/usr/bin/perl


open(OUT, ">Selected_SNP_annotation_by_chr.txt");
@selected_vcfs = `ls chr*/chr*.selected_snps.vcf`;

foreach $vcf(@selected_vcfs)
{
print "$vcf";
open(IN, "$vcf");
$version = <VCF>; $head = <VCF>;
while(<IN>)
{
(@data) = split(/\t/, $_);
if($data[7] =~ /ANN=[^\|]+\|([^\|]+)\|/){$annotation = $1; 
$chr = $data[0]; $pos = $data[1];
#print "$data[0]\t$data[1]\t$annotation\n";
$lookup{$chr}{$annotation}++; 
$annotations{$annotation}++;
}

}
close IN;

}
print OUT "Annotation";
foreach $chr(sort keys %lookup){print OUT "\t$chr";}

foreach $annotation(keys %annotations)
{
print OUT "\n$annotation";
foreach $chr(sort keys %lookup){print OUT "\t$lookup{$chr}{$annotation}";}
}

close OUT;

open(OUT, ">stats_by_bin.txt");
print OUT "Chr\tGenome\tBin\tResolved\tSNPs\n";

@outputs = `ls chr*/binned_output/out*`;
foreach $file(@outputs)
{
if($file =~ /chr([1-7])([ABD])\/binned_output\/out_bin_(\d+).txt/){$bin = $3; $chr = $1; $genome = $2;}
$last = `tail -1 $file`;
chomp $last;
($power, $resolved, $name, $data) =split(/\t/, $last);
$snps = `grep -c "^" $file`;
chomp $snps;
print OUT "$chr\t$genome\t$bin\t$resolved\t$snps\n";
}

close OUT;

open(OUT, ">number-of-snps-available-per-bin.txt");


@genotype_files = `ls chr*/chr*.vcf.single_copy.genotypes`;
foreach $file(@genotype_files)
{
chomp $file;
if($file =~ /(chr[1-7][ABD])/){$chr = $1;}
$tail = `tail -1 $chr/$chr.vcf.filtered.genotypes`;
if($tail =~ /chr[1-7][ABD]_(\d+)\t/){$length = $1;} 
open(IN, "$file");
$head = <IN>;
while($line = <IN>)
  { 
  chomp $line;
  if($line =~ /chr[1-7][ABD]_(\d+)\t/)
    {$pos = $1; $bin = int($pos/1500000); $bin_density{$chr}{$bin}++;}
  }
close IN;
}

foreach $chr(sort keys %bin_density)
{
print OUT "$chr";
$ref = $bin_density{$chr};
%hash = %$ref;
foreach $bin(sort {$a<=>$b} keys %hash){print OUT "\t$bin_density{$chr}{$bin}";}
print OUT "\n";
}
