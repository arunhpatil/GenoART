#!/usr/bin/perl
$input = $ARGV[0];
#$newinput = $ARGV[0];
$output = $ARGV[1];
$newoutput = "pre_"."$output";
#$input = "GSSPs_"."$newinput";
$gt = $output;
$gt =~ s/\.txt//;
$gt =~ s/\.csv//;
$gtfoutput = "GenoArt_"."$gt".".gtf";
open (FILE, "$input") or die $!;
open (OUT, ">$output") or die $!;
while (<FILE>)
{
	chomp ($_);
	$_ =~ s/\"//g;
	$_ =~ s/\r//g;
	$_ =~ s/\n//g;
	@lines = split (/\t/, $_);
	if ($lines[0] =~ /[P|p]eptide/ || $lines[0] =~ /[A|a]nnotated/ || $lines[0] =~ /[S|s]equence/)
	{
	}
	else
	{
		$pep = $lines[0];
		$pep =~ tr/[a-z]/[A-Z]/;
		if ($pep =~ /\./){@dots = split (/\./, $lines[0]); $pep = $dots[1]; } 
			$command = 'mysql -u arun_patil -phbp@123 -A GenoArt -e "SELECT * FROM ihp_proteogenomics WHERE peptide=\''.$pep.'\'" 2>/dev/null | grep -v "gene_info" >>'.$newoutput.'';
			system ($command);
	}
}
open (FOFI, "$newoutput") or die $!;
#open (WN, ">WhyNotMap.txt") or die $!;
while ($f = <FOFI>)
{
	chomp ($f);
	$f =~ s/\r//g;
	@data = split (/\t/, $f);
	#if ($data[9] ne "")
	#{
		$peptide_seq{$data[9]} .= "$f*";
		$data[12] =~ s/\,/\-/g;
		$ucsc_var = "<a href=\"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=".$data[0]."%3A".$data[3]."-".$data[4]."\" target=\"_blank\">$data[9]</a>";
		if ($data[11] =~ /^ENST/)
		{
			#http://www.ensembl.org/Human/Search/Results?q=ENST00000418589.1;site=ensembl;facet_species=Human
			$gene_symbol =  "<a href=http://www.ensembl.org/Human/Search/Results?q=\"$data[11]\;\"site=ensembl;facet_species=Human target=\"_blank\">$data[10]</a>";
		}
		else 
		{
			$gene_symbol =  "<a href=\"https://www.ncbi.nlm.nih.gov/gene/$data[11]\" target=\"_blank\">$data[10]</a>";
		}
		$data[13] =~ s/\,/\;/g;
		$line = "$data[1],$ucsc_var,$gene_symbol,$data[12],$data[6],$data[13],$data[14],$data[15]" ;
		#$line = "$data[1],$ucsc_var,$data[12],$data[13],$data[6],$data[9],$data[10],$data[11]" ;
#		$line = join (",", @data);
		$line_cat .= ","."$line";
	#}
	#print "$line\n";
}
$line_cat=~s/^\,//;
print "$line_cat\n";
open (FILE, "$input") or die $!;
open (GTF, ">$gtfoutput") or die $!;
while (<FILE>)
{
        chomp ($_);
        $_ =~ s/\"//g;
        $_ =~ s/\r//g;
        $_ =~ s/\n//g;
        @lines = split (/\t/, $_);
        if ($lines[0] =~ /[P|p]eptide/ || $lines[0] =~ /[A|a]nnotated/ || $lines[0] =~ /[S|s]equence/)
        {
        }
        else
        {
                $pep = $lines[0];
                $pep =~ tr/[a-z]/[A-Z]/;
                if ($pep =~ /\./){@dots = split (/\./, $lines[0]); $pep = $dots[1]; }
		if (exists $peptide_seq{$pep})
		{
			@stars = split (/\*/, $peptide_seq{$pep});
			foreach $star (@stars)
			{
				print OUT "$_\t$star\n";
				@gtfs = split (/\t/, $star);
				print GTF "$gtfs[0]\t$gtfs[1]\t$gtfs[2]\t$gtfs[3]\t$gtfs[4]\t$gtfs[5]\t$gtfs[6]\t$gtfs[7]\t$gtfs[8]\n";
			}
		}
	}
}
