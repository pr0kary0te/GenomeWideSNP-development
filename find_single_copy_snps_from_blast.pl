#!/usr/bin/perl

$max_percent = 90;
$max_length = 33;
$name = $ARGV[0];
chomp $name;
open(IN, "$name/$name.ncbi_blast.out");
open(GEN, "$name/$name.vcf.single_copy.genotypes");
open(OUT, ">$name/$name.vcf.blast.single_copy.genotypes");

while(<IN>)
{
chomp;
($seqid, $chr, $percent, $length, @other) = split(/\t/);
if($percent >$max_percent && $length > $max_length){$lookup{$seqid}++; }

}


$head = <GEN>;
print OUT $head;
while(<GEN>)
{
($id, @other) = split(/\t/);
if($lookup{$id} ==1){print OUT;}
}
