#!/usr/bin/perl


$chr = $ARGV[0];
chomp $chr;

$left = 20;
$right = 20;
$length = $left+$right+1;

$file = "$chr/$chr.vcf.single_copy.genotypes";
print "Opening $file\n";
open(FASTA, ">$chr/${chr}_single_copy_variants.fa");

$sequence ="";

$seqname = $chr; $seqname =~ s/chr//;

#Check the path points to your IWGSC chromosomes
open(CHR, "/data2/gary/IWGSCR_chromosomes_fasta/$seqname.fa");

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
    $substr= substr($sequence, ($pos - $left -1), ($length));
    $index = index($sequence, $substr);
    
    #Do some basic checks
    if($substr =~ /[CATG]{10,20}A{4,}|T{4,}|G{4,}|C{4,}[CATG]{10,20}/){$bad++;}
    else{
    print FASTA ">$name\n$substr\n"; $good++;} 
    }
  
close IN;

$total = $bad+$good;
print "$bad of $total sequences rejected due to central homoplymers\n";

