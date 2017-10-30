#!/usr/bin/perl
$input = $ARGV[0];
$output = $ARGV[1];
$not_mapped_output = "ProteinCoding_"."$output";
$mapped_output = "GSSPs_"."$output";
$fasta = $ARGV[2];
if ($input eq "" || $output eq "" || $fasta eq "")
{
	print "USAGE: perl fetch_GSSPs_IHP.pl <input_file> <output_file> <RefSeq.fasta>\n";
	exit;
}
open (FILE, "$input") or die $!;
#open (OUT, ">$output") or die $!;
open (FASTA, "./Databases/$fasta") or die $!;
while ($xf = <FASTA>)
{        
        chomp ($xf); 
        $xf =~ s/\r//g; 
        if ($xf !~ /^>/) 
        { 
                $seq .= "$xf"; 
                $seqIL .= "$xf"; 
        } 
} 
$seqIL =~ s/I/#/g; 
$seqIL =~ s/L/#/g; 
open (NOTMAPED, ">$not_mapped_output") or die $!;
open (MAPPED, ">$mapped_output") or die $!;
while (<FILE>)
{
	chomp ($_);
	$_ =~ s/\"//g;
	$_ =~ s/\r//g;
	$_ =~ s/\n//g;
	@lines = split (/\t/, $_);
	if ($lines[0] =~ /[P|p]eptide/ || $lines[0] =~ /[A|a]nnotated/ || $lines[0] =~ /[S|s]equence/)
	{
		print MAPPED "$_\n";
	}
	else
	{
		$pep = $lines[0];
		$pep =~ tr/[a-z]/[A-Z]/;
		if ($pep =~ /\./){@dots = split (/\./, $lines[0]); $pep = $dots[1]; } 
		$ilu = $pep; 
		$ilu =~ s/I/#/g; 
	        $ilu =~ s/L/#/g; 
	        #if ($seqIL !~ /$ilu/) 
	        if ($seq !~ /$pep/ && $seqIL !~ /$ilu/) 
		{
			print MAPPED "$_\n";
		}
		else
		{
			print NOTMAPED "$_\n";
		}
	}
}
