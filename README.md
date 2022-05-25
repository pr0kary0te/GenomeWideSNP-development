# GenomeWideSNP-development

This repositary contains the PERL code used to process SNPs in vcf format from ~1000 wheat skim sequence lines to generate an optimal set of genome-wide SNP markers.  The motivation was to produce a novel Axiom genotyping array.

The pipeline is mostly a set of PERL scripts, but we also use NCBI <b>BLAST+</b> (we used v2.6.0+) and the Burrows Wheeler aligner (we used <b>bwa</b> Version: 0.7.12-r1039). These need to be installed in your path, or alternatively update the scripts with an absoulte path to these tools. 

The pipline is initiated by running <b>./find_snps.pl</b>, but  file you will need to edit this file first to update some key data paths to match your local environment before you will be ready to run this pipeline. 

The first is the <b>$vcfpath</b> variable with the path to the directory containing your vcfs.  These should be named in the format chr1A.vcf, chr1B.vcf etc through to chr7D.vcf.

The second is the <b>$iwgscpath</b> variable which contains the name of the directory containing the IWGSC1.0 chromosomes in fasta format, named 1A.fa, 1B.fa etc through to 7D.fa

You can supply an optional list of files containing existing SNPs which must be included whether or not they pass filters. The files should in the present working directory and be in this comma separated format:

Chromosome (IWGSC v1),Position (IWGSC v1),35K_SNPId,Sequence
chr1A,1174865,AX-94493709,AGGAGTGATGTATCTGAAAATCTCGAGGGAGCTGA[A/G]ACCCAAGCTTCACAAGAACCTCCTGAGTTTCTGAA
chr1A,1236448,AX-94772289,GCTCGGTCGTCTTCGGCCTTGCGACCAAGCATGAT[C/T]TGGTGGATGGCAGGATGAAGAGGTGGTACCTCAGG
chr1A,1340329,AX-95211874,TCCTTCCCTGTTAATTTACCGCATGTAAGCAACAC[G/T]GGCCCACCACCTCCTCCTAAC

The files are then listed in the @selectedSNPfiles array, which you can edit as appropriate.

Finally, we added a list of core watkins landraces in our pipeline to include in the <i>core_watkins.txt</i> text file, which must be in the present working directory. There hshould be one watkins line listed on each line:
WATDE0001
WATDE0002
etc..



The find_snps.pl pipeline will run the following steps:

<b>sub_sample_elites_and_selected_landraces_only.pl</b>
This makes a copy of the original vcfs with call only for named elite lines (those not from the Watkins collection - named WAT*) along with core Watkins lines listed in (core_watkins.txt)</i>. The logic here was to balance the number of modern cultivars used for SNP selection against the much higher number of landraces available, to prevent ascertainment bias. Including some landraces is good as these will contain rare alleles which breeders will wish to bring into cultivars over time, but there is a balance to be struck with an array that will work well with existing cultivars. Note that the script refers so a file called "<i>high_heterozygosity_lines.txt</i>".  You don't need to supply this: it will be generated automatically the first time you run the pipeline.  Once you have run the pipleine, re-run it and any line identified as high heterozygosity in the first pass will be excluded in the second (they will be listed in the file so you will know which ones they were).

<b>filter_snps_from_vcf.pl</b>
SNP rows are selected only where they meet the criteria specified in the script:
$max_het_prop = 0.005; (maximum proportion of varieties with heterozygous calls)
$min_maf = 0.01; (minumum Minor allele Frequency)
$min_call_rate = 0.95; (minimum proportion of varieties with a valid hom or het call.


<b>convert_vcf_to_genotypes.pl</b>
Converts the vcf formatted SNPs passing the previous step into 0, 1, 2 coded calls for downstream analysis


<b>create_masked_chromosome_sequence.pl</b>
Prepares a copy of the IWGSC fasta format chromosomes where all SNP positions that have passed filter up to this point are masked out with Ns 20 nucleotides up and downstream. This allows us to efficiently check for SNPs with multiple off-target locations using BWA mapping using the <b>filter_out_high_copy_snps.pl<b> script.
  
generate_fasta_of_single_copy_filtered_variants.pl 
This step makes a FASTA formatted file contianing every SNP flanking sequence passing filters up to this point.  This file is then BLASTed against the local IWGSC1.0 genome and any SNPs with multiple hits filtered out with <b>find_single_copy_snps_from_blast.pl</b> which also removes SNPs with low complexity flanking regions. 

<b>get_best_snps_in_bin.pl</b> breaks the genome down into 1.5 Mb bins (by default) and uses the <b>select_minumal_markers.pl</b> script to identify the set of up to 6 SNPs (by default) that best discriminate all varieties within that bin. 
  
<b>create_affy_design_file.pl</b> gets the flanking sequence for SNPs selected in the previous step so they are ready to submit for array probe design.

  <b>check_heterozygosity_rate.pl</b> runs a post-run check on the heterozygosity rate of each variety based on the SNPs chosen by the pipleine. Any varieties with high heterozygoisty are written to the file xxx
annotate_results.pl




The next file to edit will be <b>filter_snps_from_vcf.pl</b>. In this file, we supplied three files with SNPs which were to be forced to pass the initial filters if they were found to overlap with any SNPs in the skim sequenced vcfs (we know they are good SNPs if they are present on existing arrays after all). This is important because we are definitely going to include these legacy SNPs on the new array, so we want to make sure that the newly dsigned ones are designed to work around these as well as possible - to removed redundancy. For us, these files were: <i>dart-snps.csv, cerealsdb_SNPs_from_ensembl.txt</i> and <i>breeders-additional-markers.csv</i>. You could combine these into one or add more files: just update the file list accordingly.  This step also optionally adds snpEff annotation to SNPs and this will need to be available and its path correctly specified on line 122. 

<b>create_masked_chromosome_sequence.pl</b> to provide it with a path to a local copy of the IWGSC v1.0 chromosome sequences in FASTA format. These shoud appear as $path/1A.fa, $path/1B.fa etc.  Edit the path in the $path variable and rename the chromosomes if necessary. Note that the $path directory needs to be user writable.  The purpose of this script is to make a copy of each chromosome with the canonical position of every filtered SNP masked out with N's in its flanking region. Any SNPs which then map to these masked chromosomes using bwa must have a highly similar flanking sequence in more than one chromosome location, and they will be filtered out as high copy.  


