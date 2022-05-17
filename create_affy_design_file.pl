#!/usr/bin/perl


$chr = $ARGV[0];
chomp $chr;

$leftflank = 100;
$rightflank = 100;

$existing = "$chr/affy_designs_for_$chr.txt";
if(-e $existing){`cp $chr/affy_designs_for_$chr.txt $chr/affy_designs_for_$chr.txt.bak`;}

open(DESIGN, ">$chr/affy_designs_for_$chr.txt");
print DESIGN "Name\tProbeseq\tMatches\tRank\n";

 %used =();
 @files = `ls $chr/binned_output/out*`;

foreach $file(@files)
  {
  @lines = `cat $file`;
  $i = 0;
  if($file =~ /out_bin_(\d+)/){$bin = $1;}
   foreach $line(@lines)
    {
    $i++;
    ($score, $resolved, $name, $haplotype) = split(/\t/, $line);
    $used{$name}++;
    $rank{$name} = $i;
    $name2bin{$name} = $bin;
    }  
   }

open(HEAD, "$chr/$chr.vcf.header");
while(<HEAD>){$header .= $_;} 

$sequence ="";

$seqname = $chr; $seqname =~ s/chr//;

open(CHR, "/data2/gary/IWGSCR_chromosomes_fasta/$seqname.fa");

$head = <CHR>;
#print "$head";

while(<CHR>)
{
chomp; 
$sequence.=$_;
}
close CHR;



open(VCF, ">$chr/$chr.selected_snps.vcf");
print VCF "$header";
open(IN, "/home/bzglab/rdsf/elite_plus_watseq_raw_vcf/$chr.vcf");
while(<IN>)
  {
  # chr1A	1569	chr1A_1569	G	A	 
  ($chr, $pos, $name, $ref, $alt, @other) = split(/\t/, $_);
  if($used{$name} > 0)
    {
    print VCF;
    $snp = substr($sequence, ($pos-1),1);
    $left= substr($sequence, ($pos - $leftflank -1), ($leftflank));
    $right = substr($sequence, $pos, $rightflank);
    if($left =~ /[A-Z]+([A-Z]{20})$/){$whole = "$1$snp";}
    if($right =~ /^([A-Z]{20}).*/){$whole .= $1;} 
    $match = 0;
    $len = length($whole);
    while($sequence =~ /$whole/g){$match++;}
    $rank = $rank{$name};
    $bin = $name2bin{$name};
    print DESIGN "$name\t$left\[$ref\/$alt\]$right\t$match\t${bin}.$rank\n";
    #print "$len $whole\n$name\t$snp $left [$ref\/$alt] $right $match hits\n";
    }
  }
close IN;
close VCF;
close DESIGN;



