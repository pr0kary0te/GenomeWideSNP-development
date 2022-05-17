#!/usr/bin/perl



#Get the input file name from the command line - first argument:  e.g. ./select_minimal_markers.pl All_marker_data.txt
$infile = $ARGV[0];
chomp $infile;
#print "$infile\n";
if($infile =~ /(chr[1-7][ABD])/){$current_chr = $1;}


open(REQ, "selected_markers_from_dart_and_35K_array.csv");
$head = <REQ>;
while(<REQ>)
{
($chr, $pos, @other) = split(/\,/, $_);
if($current_chr =~ /$chr/)
 {
 $id = "${chr}_$pos";
 $required{$id}++;
 }
}
close REQ;

$n = keys %required;
print "Loaded $n required markers for $current_chr\n";



@affy_files = `ls $current_chr |grep "affy_designs_for"`;

#Need to take into account any bad SNPs blacklisted from all previous iterations, otherwise they will come back! 
foreach $affy(@affy_files)
{
chomp $affy;
#print "Reading $current_chr$affy for required SNPs\n";
open(AFFY, "$current_chr/$affy");
$head = <AFFY>;
#Set repeat character length to exclude from positions adjacent to probe SNP
$r = 4;
while(<AFFY>)
{
chomp;
$poly = 0;
($id, $probe, $count, $rank) = split(/\t/, $_);
#print "ID $id\tProbe $probe\tCount $count\tRank $rank\n";
#Make some basic checks on the probe sequence.
if($probe =~ /(A{$r}|C{$r}|T{$r}|G{$r})\[/ || $probe =~ /\](A{$r}|C{$r}|G{$r}|T{$r})/ || $probe =~ /N.{0,20}\]/ || $probe =~ /\].{0,20}N/)
   {$poly++;}
if($count >1 || $poly >0){$blacklisted{$id}++; } else{$required{$id}++;}
}
close AFFY;
}

$bn = keys %blacklisted;
$rn = keys %required;
#print "$bn blacklisted IDs and $rn required IDs loaded\n";




#The proportion of variety pairs resolved when the program will stop searching
$max_iterations = 100;
if($ARGV[3] =~/\d/){$max_iterations = $ARGV[3];}
$resolution_cutoff = 100;
if($ARGV[4] =~ /\d/){$resolution_cutoff = $ARGV[4];}

$het_weight = 2;

$maxmarkers =50000;
#it was set at 50000
#This is currently set to more than the number of markers in the input file (~ 35000) but if it is set lower then markers are prioritised.
#E.g. if set to 5000, then the top 5000 markers by MAF score will be used and the rest ignored.


# ./select_minimal_markers.pl $ARGV[0] $ARGV[1] $ARGV[2] $ARGV[3]
# ./select_minimal_markers.pl data.txt 0.001 0.3 95
$min_maf = $ARGV[1];
#it was set at 0.004
#MAF is minor allele frequency.  This is set to a low level to include as many markers as possible but exclude rare error calls.
#It needs optimising.


$min_call_rate = $ARGV[2];

#it was set at 0.6
#Ignore markers with less than this proportion of good (0, 1 of 2 ) calls.  Probably needs optimising.





#Initiate some hashes.
my %matrix = ();
my %testmatrix =();


#Open the input data file handle - the input file having been specified on the command line
open(IN, "$infile");
$original = $infile;
$infile =~ s/bin_/out_bin_/;

#print "$infile\n";
open(OUT, ">$infile");


#Parse the file header before going through all the data lines
#The file format for the header is "Marker name -> Variety1 name-> Variety2 name-> Variety3 name...."
$head = <IN>;
($id, @header) = split(/\t/, $head);
$hlen = @header;


#Start reading the data here
#The file format for the data is "Marker name -> Variety1 score  -> Variety2 score-> Variety3 score...."
while(<IN>)
{
chomp;
($id, @data) = split(/\t/, $_);
%alleles = ();
if($blacklisted{$id} <1)
{
foreach $cell(@data)
  {
  #Replace any cells which aren't 0, 1 or 2 with "x"
  if($cell !~ /^[012]$/) {$cell = "x";} 
  
  #Make a hash list of alleles observed in this current row  which aren't bad "x" calls
  if($cell !~ /x/){$alleles{$cell}++;}
  }

$thislen = @data;
#Check that the header and each data row have the same number of cells, in case of file corruption. Die if not.
if($hlen != $thislen){print "$id has length of $thislen which doesn't match header ($hlen)\n"; die;}
#print "$hlen\n";
#print "$thislen\n";

$data = join("", @data);

$n_alleles = keys %alleles;

$fails = 0;

#Count the failed calls so we can work out the good call rate
while($data =~ /x/g){$fails++;}
$callrate = ($thislen - $fails)/$thislen;


#Check the call rate is above threshold and we have more than one allele observed (i.e. polymorphic) 
# and add qualifying IDs to the %pattern2id hash. Note that if any SNPs generate identical call patterns 
# then previous IDs will get overwritten as we work through the rows: this is intentional to prevent time 
# wasting with testing SNPs with identical call patterns.


if($callrate > $min_call_rate && $n_alleles >1)
  { 
   #Check we're not overwriting a required marker with one having the same pattern..
   $prev = $pattern2id{$data};
   if($required{$prev} <1) {$pattern2id{$data} = $id;}
  }
}
}
$n = keys %pattern2id;

print "Loaded marker data for $n distinct patterns\n";
close IN;
#All done with reading data from the SNP input file now.

#Now loop over the distinct SNP patterns to organise them by MAF score. This part is only really relevant if
#$maxmarkers is set to fewer than the actual number of input markers, othewise all get used anyway regardless of MAF ordering.
 
foreach $pattern(keys %pattern2id)
 {
 #Initialise this hash each iteration time to empty it.
 %charcount = ();
 $id = $pattern2id{$pattern};
 while($pattern =~ /([012x])/g) {$charcount{$1}++;}
 $zero = $charcount{0}; 
 $one = $charcount{1};
 $two = $charcount{2};
 $count = $zero+$one+$two;
 #Logic steps to work out which is the second most common call, which we'll define as the minor allele.  
 if($one >= $zero && $zero >= $two){$minor = $zero;} 
 elsif ($zero >= $one && $one >= $two){$minor = $one;}
 elsif ($zero >= $two && $two >= $one){$minor = $two;}
  
 $maf = $minor/$hlen;
#print "$id maf is $maf: 0 $zero 1 $one 2 $two\n";
  if ($maf > $min_maf)
  {
  push @{$orderbymaf{$maf}}, $pattern; 
  $pattern2maf{$pattern} = $maf;
  }
}
print "Sorting by minor allele frequency..\n";

$count = 0;
foreach $mafscore(sort {$b<=>$a} keys %orderbymaf)
{
$ref = $orderbymaf{$mafscore};
@ary = @$ref;
$l = @ary;
if($count < $maxmarkers)
 {
 foreach $pattern(@ary){$id = $pattern2id{$pattern}; if($count < $maxmarkers){$selected{$pattern} = $id;} $count++; }
 }
}



$n = keys %selected;
print "$n distinct SNP patterns\n";


#set up n x n matrix for data
my $l = @data;

$x = 0; $y = 0;

while($x < $l){$y = 0; while($y < $l){$matrix{$x}{$y} = 0; 
#print "Matrix $l build $x $y\n"; 
$y++;} $x++;
}

print "Built initial $l x $l scoring matrix\n";
$matrix_size = ($l * ($l -1))/2;
$residual = $matrix_size;

$currentscore = 1;
while($currentscore > 0)
{
$iterations++;
%testmatrix =();
$bestscore = 0;
$first = 0;
$iteration++; print "\nIteration $iteration ";
foreach $pattern1(sort keys %selected)
   {
   $first++;
   $score = 0;
   $id = $selected{$pattern1};
   %testmatrix =();
   @pattern1 = split(//, $pattern1);
   #print "$id\t$pattern1\n";
   $i = 0; $j = 0;
   $len = @pattern1;
   while($i < $len)
      {
      $ichar = $pattern1[$i];
      $j = $i +1;
      while($j <$len)
         {
         $jchar = $pattern1[$j];
         if($matrix{$i}{$j} ==0 && $ichar ne "x" && $jchar ne "x" && $ichar ne $jchar) 
             {
             $testmatrix{$i}{$j}++; $score++;  
             }  
         $j++;
         }
      $i++;
      }
#Added the line below so that only one marker at a time can be marked as required and get up-graded score.
#if($bestscore < 9999999999999999)
#{
# if ($required{$id}>0)
# {
# $score=99999999999999999999; 
# delete $required{$id};
# print "Found required id $id and set it's score to max value.\n";
# }
#}

#New section added to weight against het calls by downgrading the marker score according to the number of hets.
$hetcount = 0;
$total_length = 0;
while($pattern1 =~ /[012]/g){$total_length++;}
while($pattern1 =~ /1/g){$hetcount++;}
$newscore = int($score * (($total_length - ($hetcount * $het_weight)) / $total_length));
#Stop search loop crashing too early if a marker with lots of  hets is encountered late in the run: only apply weighting if it doesn't result in a zero score:
if($newscore >1){$score = $newscore;}

  if($score >$bestscore)
      {
      $bestscore = $score; 
      %bestmatrix = %testmatrix; 
      $bestpattern = $pattern1;     
      }
if($required{$id} >0){  $bestscore = $score;  %bestmatrix = %testmatrix;  $bestpattern = $pattern1; delete $required{$id}; last;}  
 }


       $maf = $pattern2maf{$bestpattern}; $id = $pattern2id{$bestpattern}; 
      if($bestscore >0)
         {
         $residual -= $bestscore;
         $resolved_percent = 100* ($matrix_size - $residual)/$matrix_size;
         
         #print "$bestscore\t$resolved_percent\t$id\t$maf\t$bestpattern\n"; 
         print OUT "$bestscore\t$resolved_percent\t$id\t$bestpattern\n"; 
         }
        $currentscore = $bestscore;
        if($iterations == $max_iterations || $resolved_percent >= $resolution_cutoff){$currentscore = 0;}
        $i = 0; $j = 0;
         @pattern1 = split(//, $bestpattern);
       $len = @pattern1;
       $improved =0;
        
        while($i < $len)
         {
         $j = $i +1;
         while($j <$len)
            {
            $matrix{$i}{$j} += $bestmatrix{$i}{$j};  
            $improved += $bestmatrix{$i}{$j};
            $j++;
            }
         $i++;
         }
$bestscore =0;
%bestmatrix =();
%sortbyscore =();
}

#Cleanup by deleting the input file;
`rm $original`;



sub compare
{
my $score = 0;
my $c = $l;
my $x = 0;
my $y = 0;
while($x < $c)
{
$y = 0;
while($y < $c){if($matrix{$x}{$y} ==0 && $testmatrix{$x}{$y} >0){$score++;}$y++; }
$x++;
}
return $score;
}





