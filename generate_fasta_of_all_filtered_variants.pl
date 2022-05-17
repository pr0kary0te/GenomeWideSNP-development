#!/usr/bin/perl


$chr = $ARGV[0];
chomp $chr;

$left = 20;
$right = 20;
$length = $left+$right+1;

$file = "$chr/$chr.vcf.filtered.genotypes";
print "Opening $file\n";
open(FASTA, ">$chr/${chr}_filtered_variants.fa");

$sequence ="";

$seqname = $chr; $seqname =~ s/chr//;

open(CHR, "../IWGSCR_chromosomes_fasta/$seqname.fa");

$head = <CHR>;
#print "$head";

while(<CHR>)
{
chomp; 
$sequence.=$_;
}
close CHR;

open(IN, "$file");
$head = <IN>;
while(<IN>)
    {
    ($name, @other) = split(/\t/, $_);
    if($name =~ /chr[1-7][ABD]_(\d+)/){$pos = $1;}else{die "Bad position in $name\n";} 
    $substr= substr($sequence, ($pos - $leftflank -1), ($length));
    print FASTA ">$name\n$substr\n";
    }
  
close IN;



