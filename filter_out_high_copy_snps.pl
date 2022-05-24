#!/usr/bin/perl

$path = "/data2/gary/IWGSCR_chromosomes_fasta";
$name = $ARGV[0];
chomp $name;

if($name =~ /chr([1-7])([ABD])/){$chr = $1; $genome = $2; $name = "chr$chr$genome";}


#Generate a fasta format file of the flanking sequence around each SNP in chrXX/chrXX.vcf.filtered.genoypes

$filename = "$name/${name}_filtered_variants.fa";
if(-e $filename){} else{
`./generate_fasta_of_all_filtered_variants.pl $name`;
}


#Map reads to the masked chromosome where the canonical SNP position is masked out with Ns. 
`bwa mem $path/$chr$genome.masked.fa $name/${name}_filtered_variants.fa -t 32 >$name/${name}.bwa_mapped_to_$name.sam`;

#Get lines with no matches at all if using a bwa index file which has been masked with Ns over the canonical SNP positions
@lines = `grep AS:i:0 $name/${name}.bwa_mapped_to_$name.sam`;


#Example line:
#chr1A_36563	16	1A	275663	0	41M	*	0	0	TAGTCGGCCGCAGTACAACGGGGGATCTACCGGCAGACACG	*	NM:i:0	MD:Z:41	AS:i:41	XS:i:41
foreach $line(@lines)
  {
  @data = split(/\t/, $line);
  if($data[0] =~ /chr/){$id = $data[0]; $single_copy{$id}++;}
  }

`rm $name/${name}.bwa_mapped_to_$name.sam`;


$n = keys %single_copy;
print "$n single copy SNPs selected for $name\n";


#`rm $name/${name}.bwa_mapped_to_not_$name.sam`;

open(VCFIN, "$name/$name.vcf.filtered.genotypes");


#Create a file of SNPs with no additional mappings (single copy).
open(VCFOUT, ">$name/$name.vcf.single_copy.genotypes");
$head = <VCFIN>;
print VCFOUT "$head";

while(<VCFIN>)
{
($id, @other) =split(/\t/, $_);
if($single_copy{$id}>0){print VCFOUT;}
}


