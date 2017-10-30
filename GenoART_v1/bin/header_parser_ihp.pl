#!/usr/bin/perl
open (FILE, "header_file_ihp.txt") or die $!;
while (<FILE>)
{
	chomp($_);
	@lines = split (/\t/, $_);
	$line = join ("\t", @lines);
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	$newline = "$line\t" . "Chromosome\tSource_ID\tFeature\tStart\tEnd\tScore\tStrand\tFrame\tAttribute\tPeptide_sequence\tGene\tGene_ID\tDescription\tSplicing\tCategory\tGene_type\n";
}
open (OUT, ">header_file_temp_ihp.txt") or die $!;
print OUT $newline;
