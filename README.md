# GenomeWideSNP-development

This repositary contains the PERL code used to process SNPs in vcf format from ~1000 wheat skim sequence lines to generate an optimal set of genome-wide SNP markers.  The motivation was to produce a novel Axiom genotyping array.

The pipeline is mostly a set of PERL scripts, but we also use NCBI <b>BLAST+</b> (we used v2.6.0+) and the Burrows Wheeler aligner (we used <b>bwa</b> Version: 0.7.12-r1039). These need to be installed in your path, or alternatively update the scripts with an absoulte path to these tools. 

The pipline is initiated by running <b>./find_snps.pl</b>, but  file you will need to edit this file first to update some key data paths to match your local environment before you will be ready to run this pipeline. 

The first is the <b>$vcfpath</b> variable with the path to the directory containing your vcfs.  These should be named in the format chr1A.vcf, chr1B.vcf etc through to chr7D.vcf.

The second is the <b>$iwgscpath</b> variable which contains the name of the directory containing the IWGSC1.0 chromosomes in fasta format, named 1A.fa, 1B.fa etc through to 7D.fa

We added a list of core watkins landraces in our pipeline to include <i>(core_watkins.txt)</i>. THe logic here was to balance the number of modern cultivars used for SNP selection against the much higher number of landraces available, to prevent ascertainment bias. Including some landraces is good as these will contain rare alleles which breeders will wish to bring into cultivars over time, but there is a balance to be struck with an array that will work well with existing cultivars. Note that the script refers so a file called "<i>high_heterozygosity_lines.txt</i>".  You don't need to supply this: it will be generated automatically the first time you run the pipeline.  Once you have run the pipleine, re-run it and any line identified as high heterozygosity in the first pass will be excluded in the second (they will be listed in the file so you will know which ones they were).

The next file to edit will be <b>filter_snps_from_vcf.pl</b>. In this file, we supplied three files with SNPs which were to be forced to pass the initial filters if they were found to overlap with any SNPs in the skim sequenced vcfs (we know they are good SNPs if they are present on existing arrays after all). This is important because we are definitely going to include these legacy SNPs on the new array, so we want to make sure that the newly dsigned ones are designed to work around these as well as possible - to removed redundancy. For us, these files were: <i>dart-snps.csv, cerealsdb_SNPs_from_ensembl.txt</i> and <i>breeders-additional-markers.csv</i>. You could combine these into one or add more files: just update the file list accordingly.  This step also optionally adds snpEff annotation to SNPs and this will need to be available and its path correctly specified on line 122. 

<b>create_masked_chromosome_sequence.pl</b> to provide it with a path to a local copy of the IWGSC v1.0 chromosome sequences in FASTA format. These shoud appear as $path/1A.fa, $path/1B.fa etc.  Edit the path in the $path variable and rename the chromosomes if necessary. Note that the $path directory needs to be user writable.  The purpose of this script is to make a copy of each chromosome with the canonical position of every filtered SNP masked out with N's in its flanking region. Any SNPs which then map to these masked chromosomes using bwa must have a highly similar flanking sequence in more than one chromosome location, and they will be filtered out as high copy.  


