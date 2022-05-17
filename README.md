# GenomeWideSNP-development

This repositary contains the PERL code used to process SNPs in vcf format from ~1000 wheat skim sequence lines to generate an optimal set of genome-wide SNP markers.  The motivation was to produce a novel Axiom genotyping array.

The pipline is initiated by running ./elete_wrapper.pl, but several file will need to have their file paths updated to your local environment before you will be ready to run this pipeline.  I would like to have made everything local/relative paths but the vcf files were very large compared to our available disk space, which necessitated linking out to our research data storage system.  

The first file you need to edit is sub_sample_elites_and_selected_landraces_only.pl. Here you need to edit the $path variable with the path to your vcfs.
We also added a list of core watkins landraces here to include (file = core_watkins.txt). THe logic here was to balance the number of modern cultivars used for SNP selection against the much higher number of landraces available, to prevent ascertainment bias. Including some landraces is good as these will contain rare alleles which breeders will wish to bring into cultivars over time, but there is a balance to be struck with an array that will work well with existing cultivars. Note that the script refers so a file called "high_heterozygosity_lines.txt".  You don't need to supply this: it will be generated automatically the first time you run the pipeline.  Once you have run the pipleine, re-run it and any line identified as high heterozygosity in the first pass will be excluded in the second (they will be listed in the file so you will know which ones they were).

The next file to edit will be filter_snps_from_vcf.pl. In this file, we supplied three files with SNPs which were to be forced to pass the initial filters if they were found to overlap with any SNPs in the skim sequenced vcfs (we know they are good SNPs if they are present on existing arrays after all).


