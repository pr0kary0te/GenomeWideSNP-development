# GenomeWideSNP-development

This repositary contains the PERL code used to process SNPs in vcf format from ~1000 wheat skim sequence lines to generate an optimal set of genome-wide SNP markers.  The motivation was to produce a novel Axiom genotyping array.

The pipline is initiated by running ./elete_wrapper.pl, but several file will need to have their file paths updated to your local environment before you will be ready to run this pipeline.  I would like to have made everything local/relative paths but the vcf files were very large compared to our available disk space, which necessitated linking out to our research data storage system.  

The first file you need to edit is 
