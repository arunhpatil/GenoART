#!/usr/bin/perl

$gtf_file = $ARGV[0];
@dots = split (/\./, $gtf_file);

$input = $gtf_file;
$input =~ s/.txt/_categorized.gtf/g;
print "$gtf_file\n$input\n";
#$input = "$dots[0]"."_categorized.gtf";

$output = "$dots[0]"."_nonRed.gtf";
open (GTF, "$input") or die $!;
#open (GTF,"gencode_all_annotation_reverse_No_MC_040617.gtf") or die $!;
while ($gt = <GTF>)
{
        chomp($gt);
        $gt =~ s/\r//g;
		#print "$gt\n";
        @lines = split (/\t/, $gt);
        @semi = split (/\"/, $lines[8]);
	$redudancy_check = "$lines[0]\_$lines[6]\_$semi[3]";
	$occurance{$redudancy_check}++;
	#$pep_key{$semi[3]} .= 
	#$pepseq{$semi[3]} .= "$gt*";
}
while(($x, $y) = each %occurance)
{
	@under = split (/\_/, $x);
	$pepred{$under[-1]} .= "$under[0]\_$under[1];";
	#print "$x\t$y\n";
}
while (($pra, $prb) = each %pepred)
{
	@semi = split (/\;/, $prb);
	if ($#semi == 0)
	{
		$required{$pra} = "$prb";
		#print "$pra\t$prb\n";
	}
	#else
	#{
	#	print "-->$pra\t$prb\n";
	#}
}
open (OUT, ">$output") or die $!;
open (GTF, "$input") or die $!;
while ($gt = <GTF>)
{
        chomp($gt);
        $gt =~ s/\r//g;
        @lines = split (/\t/, $gt);
        @semi = split (/\"/, $lines[8]);
	if (exists $required{$semi[3]})
	{
		print OUT "$gt\n";
	}
}
#chr22	101011	CDS	29231621	29231641	.	+	.	gene_id "ENST00000404755.7"; transcript_id "AHRPQGR";
#chr22	101012	CDS	29259298	29259396	.	+	.	gene_id "ENST00000404755.7"; transcript_id "AGGSLPPPPRPNLGEPSCAPPSLLSSAGLPTPR";
#chr22	100113	CDS	29233461	29233468	.	+	.	gene_id "ENST00000404755.7"; transcript_id "AHGSPWASWPHR";
