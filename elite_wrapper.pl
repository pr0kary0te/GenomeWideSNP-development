#!/usr/bin/perl


#Pipeline to identify the optimal set of SNP markers accross the wheat genome, starting with a set of VCF files called on skim sequencing mapped to the IWGSCv1 reference.

#File paths that need to be set up according to your local environment

#The location of your wheat vcf files. They should be in the following folder and named as: chr1A.vcf, chr1B.vcf ... chr7D.vcf
$vcfpath = "/home/bzglab/rdsf/elite_plus_watseq_raw_vcf/";

#IWGSC chromosomes path
$iwgsc = "/data2/gary/IWGSCR_chromosomes_fasta/";

#List of files containing existing SNPs which must be included whether or not they pass filters - this is optional. The files should in the pw and be in this comma separated format:

#Chromosome (IWGSC v1),Position (IWGSC v1),35K_SNPId,Sequence
#chr1A,1174865,AX-94493709,AGGAGTGATGTATCTGAAAATCTCGAGGGAGCTGA[A/G]ACCCAAGCTTCACAAGAACCTCCTGAGTTTCTGAA
#chr1A,1236448,AX-94772289,GCTCGGTCGTCTTCGGCCTTGCGACCAAGCATGAT[C/T]TGGTGGATGGCAGGATGAAGAGGTGGTACCTCAGG
#chr1A,1340329,AX-95211874,TCCTTCCCTGTTAATTTACCGCATGTAAGCAACAC[G/T]GGCCCACCACCTCCTCCTAAC


@selectedSNPfiles = ("dart-snps.csv","cerealsdb_SNPs_from_ensembl.txt","breeders-additional-markers.csv");





#This pipeline can either be given a specific chromosome to analyse as a single command line argument, e.g. 1A, otherwise it will work through all 21 (1A - 7D)
$chr = 0;
$name = $ARGV[0];
chomp $name;
if($name =~ /([1-7])([ABD])/){$chromosomes{$1}++; $genomes{$2}++; print "Working on chr$1$2\n";}


#No chromosome supplied so add all possible chromosomes from 1-7 ABD to the to do list.
else{
$genomes{A}++;
$genomes{B}++;
$genomes{D}++;

$chromosomes{1}++;
$chromosomes{2}++;
$chromosomes{3}++;
$chromosomes{4}++;
$chromosomes{5}++;
$chromosomes{6}++;
$chromosomes{7}++;
}


#Comment in or out sections to run: this will allow individual sections to be re-run with different parameters without starting from scratch.
#If running for the first time, you would expect all sections to be activated by uncommenting them.

$section{1}++;  #Start from full vcf files and get the data only for the varieties we want
$section{2}++;  #Filter SNPs in the vcf file by quality and then convert these to the genotype format
$section{3}++;  #Filter by copy number using BWA mapping to just the canonical or to all homeologous chromosomes.
$section{4}++;   #Do a BLAST search to check that no multi copy SNPs got through just because tehy were not a good enough match for BWA to pick up
$section{5}++;   #Do final SNP calling.


#Set up two foreach loops, one for chromosomes and another for genomes: this will then work through 1A, 1B, 1D then 2A, 2B, 2D and so on until 7D
foreach $chr(sort {$a<=>$b} keys %chromosomes)
{
foreach $genome (sort keys %genomes)
 {
 $name = "chr$chr$genome"; 
  $chromosome = $name;
  $chromosome =~ s/chr//;

if($section{1}>0)

   {
   print "Working on $name\nSub-sampling elites.";
@out =  `./sub_sample_elites_and_selected_landraces_only.pl $name $vcfpath`;
#Note this script checks for and includes core Watseq lines in core_watkins.txt but exludes all others.
#It then looks for bad (high het vars identified with check_heterozygosity_rate.pl in file high_heterozygosity_lines.txt
print "@out";
   }



if($section{2}>0)
{
 print "Filename $name/$name.vcf filtering by quality and pre-existing array presence: output is $name/$name.vcf.filtered\n";

$file = "$name\$name.vcf";
@out =  `./filter_snps_from_vcf.pl $file @selectedSNPfiles`;
print "@out";

`./convert_vcf_to_genotypes.pl $name/$name.vcf.filtered`;
}





if($section{3}>0)

{

#Prepare a version of the current chromosome with canonical SNP poistions and flanking regions masked with N's - to ensure any BWA mappings are to non-canonica$
  $chromosome = $name;
  $chromosome =~ s/chr//;
  #Check the path here points to a directory with the masked chromosome sequence file created above
  $bwafile = "$iwgsc/$chromosome.masked.fa";


     `./create_masked_chromosome_sequence.pl $name $iwgsc`;
      print "Creating bwa index of $bwafile\n";
      $done = "$bwafile.pac";
      if (-e $done){} else{`bwa index $bwafile`;}


  print "Running filter_out_high_copy_snps.pl $name\n";

 

  

#This script will only filter out SNPs with flanks mapping to >1 locus in the cnaonical chromosome so is less severe
$out = `./filter_out_high_copy_snps.pl $name`;
print "$out\n";

}


if($section{4}>0)
{
#Also perform a blast search (more sensetive) to take out multi-copy SNPs 
#This script also exlcudes homopolymer-contianing probes.  Input is $name/$name.vcf.single_copy.genotypes";
#Output is $name/${name}_single_copy_variants.fa");

`./generate_fasta_of_single_copy_filtered_variants.pl $name`;
#And then BLAST the resulting fasta against the corresponding IWGSC chromosome to check for the number of hits to the canonical chromosome.

#Comment out the line below to blast against only the canonical chromosome, otherwise the blast will be to an alias file containing all three homeologs 
$chromosome =~ s/[ABD]//;

$out = `blastn -query $name/${name}_single_copy_variants.fa -db ../IWGSCR_chromosomes_fasta/$chromosome.fa -num_threads 48 -outfmt 6 -evalue 1e-10 -out $name/$name.ncbi_blast.out`;




print "$out";

#And filter this output to include only single copy (by BLAST) genotypes as input for SNP finding.
`./find_single_copy_snps_from_blast.pl $name`;
`mv $name/$name.vcf.single_copy.genotypes $name/$name.vcf.single_copy.genotypes.bak`;
`mv $name/$name.vcf.blast.single_copy.genotypes $name/$name.vcf.single_copy.genotypes`;
}


if($section{5} >0)
{


$out = `./get_num_resolved_per_physical_bin_all_vars_wrapper.pl $name`;
print "$out";
#`rm $name/$name.vcf.filtered`; 
#`rm $name/$name.vcf.filtered.genotypes`;


#Do some cleanup:
`rm $name/$name.*.sam`;
`rm $name/$name.*.fa`;

system("./create_affy_design_file.pl $name &");
}

#Creates the file affy_designs_for_chrXX.txt which designs the flanking sequence for the probes and checks the number of times this
# (or a slightly shorter version) matches the genome. In the next iteration of this wrapper (if needed) any probes with >1 match 
#will be excluded by the select_minimal_markers script and any good ones forced into the design. 

}
}

$out = `./check_heterozygosity_rate.pl`;
#Creates a new list of high heterozygoisty lines in case any new ones show up and need to be removed from the next iteration of this design process.

#Finally run post_SNP_selection_summary2.pl manually, make a record of the results and then re-run this pipeline one or more times to re-design with SNPs which 
#are not high copy number and don't have homopolymers close to the central SNP etc as defined in the select_minimal_markers.pl script.


#Do some stats
`./annotate_results.pl`;
