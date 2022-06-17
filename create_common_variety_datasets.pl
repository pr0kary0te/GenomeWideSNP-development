#!/usr/bin/perl

foreach $file(@ARGV)
{
chomp $file;
$head = `head -1 $file`;
chomp $head;
($id, @data) = split(/[\t\,]/, $head);

foreach $var(@data)
 {
 $var =~ s/^[A-Z][[0-9][0-9]_\d+_//;
 $var =~ s/^[A-Z][[0-9][0-9]_//;
 $var =~ s/[^A-Za-z0-9\_]//g;

 $var =~ s/Watkin*(_[0-9]+)\./Watkins${1}_/;
 # $var =~ s/(Watkin*_[0-9]+)\./${1}_/;
 $lookup{$var}++;
 }
}

$n = keys %lookup;
$files = @ARGV;
print "$n varieties loaded accross $files files\n";
foreach $var(keys %lookup){$n = $lookup{$var};  if($n>= $files) {print "$var is present in $files files\n"; $good++;}}
print "$good varieties were shared between all $files files\n";

foreach $file(@ARGV)
{
chomp $file;
open(IN, "$file");
open(OUT, ">common$file");
$head = <IN>;
($id, @header) = split(/[\t\,]/, $head);
$datacol =0;
print OUT "ID";
foreach $var(@header)
 {
 $var =~ s/^[A-Z][[0-9][0-9]_\d+_//;
 $var =~ s/^[A-Z][[0-9][0-9]_//;
 $var =~ s/[^A-Za-z0-9\_]//g;
 $var =~ s/Watkin*(_[0-9]+)\./Watkins${1}_/;
 $count = $lookup{$var};
 $usedvar{$var}++; 
 if($count >= $files && $usedvar{$var} ==1){$common{$datacol}++; print OUT ",$var";}  
 $datacol++;
 }
print OUT "\n";

while(<IN>)
{
chomp;
($id, @data) = split(/[\t\,]/);
$datacol =0;
print OUT "$id";
foreach $cell(@data){if($common{$datacol} >0){print OUT "\,$cell";} $datacol++;}
print OUT "\n";
}
close IN;
close OUT;
}
