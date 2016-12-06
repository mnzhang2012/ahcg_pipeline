#! /usr/bin/perl -w
use strict;

open(IN,"dcm_gene_list.txt");
my $com = "grep \"";
foreach(<IN>){
    chomp;
    my @t = split("\t",$_);
    $com = "$com$t[1]\\|";
	  
}
chop $com;
chop $com;
$com = "$com\" input_files/ref_genome/hg19_refGene.txt > temp/bed_list.txt";
#print $com;
system($com);
close IN;

open(BED,"temp/bed_list.txt");
open(OUT,">temp/temp.bed");


foreach(<BED>){
	my @st;
	my @end;
	my @t;
    @t = split("\t",$_);
    @st = split(",",$t[9]);
    @end = split(",",$t[10]);
    
    my $i = 0;

	while($i < @st){
	    my $start = $st[$i]-20;
	    my $end = $end[$i]+20;
        $t[2] =~ s/chr//;
        print OUT "$t[2]\t$start\t$end\t$t[1]\t$t[8]\t$t[3]\n";
        $i++;
	}
}

close BED;
close OUT;

$com = "sort -k1n -k2 -o temp/temp.bed temp/temp.bed";
system($com);

open(BED,"temp/temp.bed");
open(OUT,">dcm_exon_list.bed");

foreach(<BED>){
    print OUT "chr";
    print OUT $_;
}

close BED;
close OUT;

$com = "rm temp/temp.bed";
system($com);

exit;
