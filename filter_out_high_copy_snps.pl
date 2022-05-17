#!/usr/bin/perl


$name = $ARGV[0];
chomp $name;

if($name =~ /chr([1-7])([ABD])/){$chr = $1; $genome = $2; $name = "chr$chr$genome";}


#Generate a fasta format file of the flanking sequence around each SNP in chrXX/chrXX.vcf.filtered.genoypes

$filename = "$name/${name}_filtered_variants.fa";
if(-e $filename){} else{
`./generate_fasta_of_all_filtered_variants.pl $name`;
}

#Map the reads in the resulting FASTA file to the other two genomes first and remove any that match
#E.g. map chr1A sequences to 1B and 1D in /data2/gary/IWGSCR_chromosomes_fasta/not_chr1A.fa


`bwa mem /data2/gary/IWGSCR_chromosomes_fasta/$chr$genome.masked.fa $name/${name}_filtered_variants.fa -t 32 >$name/${name}.bwa_mapped_to_$name.sam`;
#Cleanup fasta file as no longer needed.
#`rm $name/${name}_filtered_variants.fa`;
#Now select single mapping reads only to exlude those hitting multiple locations in the same genome as the SNP

#Get only lines with a single mapping (use with regular SAM genome mapping
#@lines = `grep -v XA:Z $name/${name}.bwa_mapped_to_$name.sam`;

#Or get lines with no matches at all if using a bwa index file which has been masked with Ns over the canonical SNP positions
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
open(VCFOUT, ">$name/$name.vcf.single_copy.genotypes");
$head = <VCFIN>;
print VCFOUT "$head";

while(<VCFIN>)
{
($id, @other) =split(/\t/, $_);
if($single_copy{$id}>0){print VCFOUT;}
}


