#!/usr/bin/perl

$minmaf = 0.01;
$bins = 20;
$file = $ARGV[0];
chomp $file;

open(IN, "$file");
open(OUT, ">final_snp_distribution_${file}min_maf-$minmaf.txt");


while(<IN>)
{
chomp;
($id, $chr, $pos, @other) = split(/[\t\,]/);
$maf = $other[$#other];
if($chr =~ /[1-7][ABD]/ && $pos =~ /^\d+$/ && $maf >= $minmaf){$matrix{$chr}{$pos}++; $used_probes++; if($pos > $max{$chr}){$max{$chr} = $pos;}}
}
close IN;
print "$used_probes exceeded min maf of $minmaf\n";
foreach $chr(sort keys %matrix)
{
$max = $max{$chr};
$ref = $matrix{$chr};
%hash = %$ref;


foreach $pos(sort {$a<=>$b} keys %hash)
{
$bin = int(($pos/$max) * $bins);
$final{$chr}{$bin}++;
}

}

$bin = 0;
while($bin < $bins)
{
$bin++;
print OUT "\t$bin";
}
print OUT "\n";
foreach $chr(sort keys %final)
{
print OUT "$chr";
$bin = 0;
while($bin < $bins)
{
print OUT "\t$final{$chr}{$bin}";
$bin++;
}
print OUT "\n";
}
