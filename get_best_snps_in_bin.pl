#!/usr/bin/perl

#Number of concurrent threads to allow
$max_threads = 48;

#Example command line:  
# ./get_best_snps_in_bin.pl chr1A/chr1A.ann.vcf

#Set the bin size within which to find the selected number of minimal markers
$bin_size = 1500000;

#Number of SNPs to find for each bin - the minimal SNP program will stop at this point and we can then see what proportion of SNPs were resolved.
$max_iterations = 6;
$resolution_cutoff = 75;


#We want to start with the filtered VCF files which have had het SNPs removed 
#using the ./filter_snps_from_vcf.pl script with these parameters:
#$max_het_prop = 0.005; (maximum proportion of het (0/1) call allowed among all varieties
$min_maf = 0.01; 
#(minor allele frequency min)

$min_call_rate = 0.95; 
#(The proportion of SNPs with data)

#Note the filtered output files have also been annoted with SNPeFF
#The output file with all SNPeFF annoted SNPs is in the name format chr6A.ann.vcf - the original unfiltered files are huge and not available here
#The subset of only non-synonymous SNPs is in the name format chr6A.non_synonymous.vcf

#We can process either the complete filtered list of just the non-synonymous list.


$name = $ARGV[0];
chomp $name;
`mkdir $name/binned_output`;

print "Input file is: $file\n";
`rm $name/binned_output/*`;


$bin = 1;
   print "Trying to open $file\n"; 
   close IN;
 
$file = "$name/$name.vcf.single_copy.genotypes";

$lines = `grep -c "^" $file`;
chomp $lines;
print "$file has $lines lines\n";
$lines--;

$first_data_line = `head -2 $file|tail -1`;
$last_data_line = `tail -1 $file`;


if($first_data_line =~ /^chr[1-7][ABDT]_(\d+)\t/){$start = $1;} else{die "Can't get start coordinate from $first_data_line\n";}
if($last_data_line =~ /^chr[1-7][ABDT]_(\d+)\t/){$end = $1;} else{die "Can't get start coordinate from $last_data_line\n";}
$range = $end - $start;


if($bins >0){$bin_size = int($range/$bins);}
print "range is $range bins $bins bin size is $bin_size\n"; 
  open(IN, "$file");

   $head = <IN>;
   
   #Make an output directory if it doesn't already exist
   close OUT;
   open(OUT, ">$name/binned_output/bin_$bin.txt");
   print OUT "$head";
   $locus = 0;
   while(<IN>)
     {
     $tally++;
     if($_ =~ /^chr[1-7][ABDT]_(\d+)\t/){$pos = $1;}
     $cumulative = ($pos - $locus); 
     if($cumulative >= $bin_size)
        {
        #$jobs++;
        $locus= $pos; $cumulative = 0; close OUT;  
        #Command line to run minimal SNP program with min MAF of 0.001, min call rate of 0.5 and stop searching for markers when 95% of var pairs are resolved
        
        #Check there are not too many unprocessed input files before launching a new thread..
 
        @files = `ls $name/binned_output/bin*`;
        $threads = @files;
        if($threads >$max_threads){while($threads >$max_threads){ @files = `ls $dir/binned_output/bin*`; $threads = @files; sleep 10;}}     
   
        system("./select_minimal_markers.pl $name/binned_output/bin_$bin.txt $min_maf $min_call_rate $max_iterations $resolution_cutoff&");
        $bin++;
        open(OUT, ">$name/binned_output/bin_$bin.txt"); print OUT "$head";
       #if($jobs ==48){$jobs = 0; sleep 300;}  
      }

     print OUT;
     }

   close OUT;

#Cleanup the last file - which may be a smaller bin
system("./select_min_markers.pl $name/binned_output/bin_$bin.txt 0.001 0.5 $max_iterations $resolution_cutoff");


#`rm $file`;
