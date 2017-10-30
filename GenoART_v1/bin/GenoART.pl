#!/usr/bin/perl
################################
### WIthout Subroutine  Extension to new5.pl###
################################



unless(scalar (@ARGV) == 4 )
{
	print "Give the proper format: perl GenoART.pl 0/1/2_MC Enzyme_name input.fasta\n";
	exit;
}
else
{
	#chomp @ARGV;
	$infile=$ARGV[2];
	#$outfile=$ARGV[1];
	$cleavage=$ARGV[0];
	$enzyme=$ARGV[1];
	$genetic_code_tt = $ARGV[3];
}

$out_s2s=$out_gtf=$in_file=$infile;
($type)=(split /\./,$infile)[-1];
$in_file=~s/\.$type/\_op\.$type/;

#$out_gtf=$outfile."_".$cleavage."MC.gtf";

$out_gtf=~s/\.$type/\_pep\.gtf/;
$out_s2s=~s/\.$type/\_pep\.fasta/;


open FILE, "$in_file" or die $!;
open (GENETIC, "bin/genetic_code/$genetic_code_tt") or die $!;


while ($genetic = <GENETIC>)
{
	chomp($genetic);
	@tab = split (/\t/, $genetic);
	$codon{$tab[0]} = "$tab[1]";
}
while (<FILE>)
{
	chomp ($_);
	if ($_ =~ /^\>/)
	{
		$x = $_;
		$x =~ s/\>//;
		@xarr = split (/\s/, $x); 
		$gene{$xarr[1]} = "$xarr[2]\|$xarr[3]"; 
	}
	else
	{
		$seq{$xarr[1]} .= "$_";  
		$seq_s2s{$xarr[1]} .= "$_"; 
	}	
}
$giaccession=0;
##########################


#$out_fasta=$outfile."_".$cleavage."MC.fasta";



#open (ALL, ">$out_fasta") or die $!;
open (GTF, ">$out_gtf") or die $!;
open (OUT, ">$out_s2s") or die $!;
#################################
while (($a,$b) = each %seq)  
{
	
	@sequence = split (//, $b);
	$seq_size = @sequence; 
	#$mycount++;
	#print "$mycount\n";
	#frame($a,\@sequence,\%codon,$giaccession,\%gene,$enzyme,$cleavage,$seq_size);
	#sub frame{
	#my ($a,$arrref1,$hashref1,$giaccession,$hashref2,$enzyme,$cleavage,$seq_size)=@_;
	#my @sequence= @{$arrref1};
	#my %codon= %{$hashref1};
	#my %gene= %{$hashref2};
	
	#my $protein="";
	#my @split_seq="";
	my @from_mc1=@seq_split_mc="";
	#print ALL "\nfirst @split_seq\n";
	#print ALL "\nProtein $protein\n";

	$seq_size=scalar(@sequence);

	for ($j =0; $j<=2; $j++)       # Creatign three frames j= 0,1,2
	{
	#print "\nframe $j\n";
	#system ("date");
		for ($i=$j; $i<= $seq_size; $i=$i+3) # take the three ltters and compare it with codon
		{
			$triplet = "$sequence[$i]"."$sequence[$i+1]"."$sequence[$i+2]"; #mergin the three letter for one codon
			$protein .= "$codon{$triplet}"; # Fetch the amino acids for the above three letter
		}
		$prt_len = length $protein; # taking the length of concatinated amino acids
		$seq = $protein; #Storing the ptoeins in seq
		$for_pep_pos = $protein;
		$count = $j+1; 
		if($enzyme eq "trypsin")
		{
			$arrref=trypsin($seq);
			@split_seq=@{$arrref};
			#print ALL "Trypsin @split_seq\n";
		}
		elsif($enzyme eq "chymotrypsin")
		{
			$arrref=chymotrypsin($seq);
			@split_seq=@{$arrref};
		}
		elsif($enzyme eq "argC")
		{
			$arrref=argC($seq);
			@split_seq=@{$arrref};
		}
		elsif($enzyme eq "lysC")
		{
			$arrref=lysC($seq);
			@split_seq=@{$arrref};
		}
		elsif($enzyme eq "lysN")
		{
			$arrref=lysN($seq);
			@split_seq=@{$arrref};
		}
		elsif($enzyme eq "aspN")
		{
			$arrref=aspN($seq);
			@split_seq=@{$arrref};
		}
		elsif($enzyme eq "gluC")
		{
			$arrref=gluC($seq);
			@split_seq=@{$arrref};
		}

		if($cleavage == 1)
		{
			$arrref1=mc1(\@split_seq);
			@seq_split_mc=@{$arrref1};
			push @seq_split_mc,@split_seq;
			@split_seq=@seq_split_mc;@seq_split_mc="";
		}
		elsif($cleavage == 2)
		{
			#print ALL "\n\nNew @split_seq\n\n";
			$ref_from_mc1=mc1(\@split_seq);
			@from_mc1=@{$ref_from_mc1};
			push @from_mc1,@split_seq;
			#print ALL "santosh4 @from_mc1\n";
			$arrref1=mc2(\@split_seq);
			@seq_split_mc=@{$arrref1};
			#print ALL "New2 @seq_split_mc\n";
			push @seq_split_mc,@from_mc1;
			@split_seq=@seq_split_mc;@seq_split_mc=@from_mc1="";

		}
                #@split_seq=@seq_split_mc;@seq_split_mc="";
		#print "@split_seq\n";
		#$giacession
		#print ""
		$uref=unique(\@split_seq);
		@split_seq=@{$uref};
		#print ALL "santosh3 @split_seq\n";
	#print "\nstart tryptic\n";
	#system ("date");		
		#tryptic(\@split_seq,$giaccession,$a,\%gene,$count,$seq_size);
		##############
		#sub tryptic
#{

		#my ($arrref1,$giaccession,$a,$hashref1,$count,$seq_size)=@_;
		#my @split_seq=@{$arrref1};
		#my %gene=%{$hashref1};
	#$giaccession++;
	#print "$giaccession\n";
	#print "$count\t@split_seq\n";
	#print "\nstart tryptic Inside\n";
	#system ("date");
	#print "@split_seq\n";		

		foreach $var (@split_seq)
		{
			if ($var !~ /\*/) # if the tryptic peptide not contains *
			{	
				$var_str =$var;
				$giaccession_new++;
				if ($var_str =~ /^M/)
				{
					for ($actyl=0; $actyl<3; $actyl++)
					{
						#@pos=();
						#print "santosh $var\n";
						if($actyl == 0)
						{
							$var =$var_str;
						}
						else
						{
							substr $var_str, 0, 1, "";
							$var = $var_str;
						}
						#$seq_len = length ($var
						#$giaccesion_new=
						#coordinate_fetcher($var,$giaccession_new,$a,\%gene,$count,$seq_size,$for_pep_pos);
						##################################################start fo coordinate fetcher
						
						#sub coordinate_fetcher
						#{
						#my ($var,$giaccession,$a,$hashref1,$count,$seq_size,$for_pep_pos)=@_;
			#print "\n<---- $count ---->\n";
						@pos=();
	#system ("date");
						#my %gene=%{$hashref1};
	#system ("date");	
						$seq_len = length ($var);
						#if ($seq_len >= 7)
						if ($seq_len >= 6 && $seq_len <= 35)
						{
							$giaccession++;
							$giaccession_1m = $giaccession;
					#print "$giaccession_1m\n";
					@coordis = split (/\|/, $gene{$a});
					$strandis = chop ($coordis[0]);
					if ($strandis eq "-")
					{
						$temp_pep_pos = reverse $for_pep_pos;
						$var = reverse $var;
						$result = index($temp_pep_pos, $var);
						$start = (($result+1)*3)-3; #((index is 0 based therefore +1)it is multiplied by 3 to get the triplet count and ) substracted by 2 because if triplet of M is ATG the start of the sequence will 3-2 i.e., ATG - TG = A;
						if ($count == 2)
						{
							$start = $start+2;
						}
						elsif ($count == 3)
						{
							$start = $start+1;
						}
					}
					else
					{
						$result = index($for_pep_pos, $var);
						$start = (($result+1)*3)-3; #((index is 0 based therefore +1)it is multiplied by 3 to get the triplet count and ) substracted by 2 because if triplet of M is ATG the start of the sequence will 3-2 i.e., ATG - TG = A;
						if ($count == 2)
						{
							$start = $start+1;
						}
						elsif ($count == 3)
						{
							$start = $start+2;
						}
					}
					$end = ($start+($seq_len*3))-1;
					$junc_end = $end - $start; $end_bckup = $end;
					@coordis = split (/\|/, $gene{$a});
					$strandis = chop ($coordis[0]);
					if ($coordis[1] !~ /\,/)
					{
						@exons_cors = split (/\-/, $coordis[1]);
						if ($strandis eq "+" || $strandis eq "\.")
						{
							$genomic_start = $exons_cors[0]+$start;
							$genomic_end = $genomic_start+($end-$start);
	#						print ALL ">gi\|10101$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$genomic_start\-$genomic_end\n";
        #	        		                print ALL "$var\n";
							print GTF "$coordis[0]\t10101$giaccession_1m\tCDS\t$genomic_start\t$genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$var\";\n";
						}
						elsif ($strandis eq "-")
						{
							$re_seq_size = $seq_size%3;
							if ($re_seq_size == 1)
							{
								if($count == 2)
								{		
									$genomic_start= ($exons_cors[0]+$start+1)-3;
									$genomic_end = $genomic_start+($end-$start);
								}
								else
								{
									$genomic_start= $exons_cors[0]+$start+1;
									$genomic_end = $genomic_start+($end-$start);
								}	
							}
							elsif ($re_seq_size == 2)
							{
								if($count == 2)
								{
									#print "santosh\n";
									$genomic_start= ($exons_cors[0]+$start+2)-3;
									$genomic_end = $genomic_start+($end-$start);
								}
								elsif($count == 3)
								{		
									$genomic_start= ($exons_cors[0]+$start+2)-3;
									$genomic_end = $genomic_start+($end-$start);
								}
								else
								{
									$genomic_start= $exons_cors[0]+$start+2;
									$genomic_end = $genomic_start+($end-$start);
								}	
							}
							else
							{
								$genomic_start= $exons_cors[0]+$start;
								$genomic_end = $genomic_start+($end-$start);
							}
							#print "behera\n";
	#						print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$genomic_start\-$genomic_end\n";
							$var_var = $var;
        		        	                $rev_var = reverse $var_var;
        #		        	                print ALL "$rev_var\n";
							print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$genomic_start\t$genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
						}
					}
					else 
					{
						@comma = split (/\,/, $coordis[1]);
						$start_bckup = $start;
						$index_comma=0;
						if ($strandis eq "+" || $strandis eq "-")
						{
							$flag=0;
							foreach $exonpair (@comma)
							{
								$index_comma++;
								@exons_cors = split (/\-/, $exonpair);
								if (($exons_cors[0]+$start) < $exons_cors[1])
								{
									$new_genomic_start = $exons_cors[0]+$start;
									$new_genomic_end = $new_genomic_start+($end-$start_bckup);
									if ($new_genomic_end <= $exons_cors[1])
									{
										if ($strandis eq "+")
										{
	#										print ALL ">gi\|10101$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
	#			        	        		                print ALL "$var\n";
											print GTF "$coordis[0]\t10101$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$var\";\n";
										}
										elsif ($strandis eq "-")
										{
											$re_seq_size = $seq_size%3;
											if ($re_seq_size == 1)
											{
												if($count == 2)
												{
													$new_genomic_start=($new_genomic_start+1)-3;
													$new_genomic_end=($new_genomic_end+1)-3;
												}
												else
												{
													$new_genomic_start=$new_genomic_start+1;
													$new_genomic_end=$new_genomic_end+1;
												}		
	#											print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
												$var_var = $var;
					        		        	                $rev_var = reverse $var_var;	
	#				        		        	                print ALL "$rev_var\n";
												print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
											}
											elsif ($re_seq_size == 2)
											{
												if($count == 2)
												{
													$new_genomic_start=($new_genomic_start+2)-3;
													$new_genomic_end=($new_genomic_end+2)-3;
												}
												elsif($count == 3)
												{
													$new_genomic_start=($new_genomic_start+2)-3;
													$new_genomic_end=($new_genomic_end+1)-2;
												}
												else
												{
													$new_genomic_start=$new_genomic_start+2;
													$new_genomic_end=$new_genomic_end+2;
												}	
												$var_var = $var;
					        		        	                $rev_var = reverse $var_var;	
	#											print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
	#				        		        	                print ALL "$rev_var\n";
												print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
											}
											else
											{
	#											print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
												$var_var = $var;
					        		        	                $rev_var = reverse $var_var;	
	#				        		        	                print ALL "$rev_var\n";
												print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
											}
										}
										last;
									}
									else
									{
										$p = $junc_end - ($exons_cors[1]-$new_genomic_start)-1;
										for ($i = $index_comma; $i <= $#comma; $i++)
										{
											@jun_exons_cors = split (/\-/, $comma[$i]);
											if (($jun_exons_cors[0]+$p) <= $jun_exons_cors[1])
											{
												#JUNCTIONS
												$new_genomic_end = $jun_exons_cors[0]+$p;
												if ($strandis eq "+")
		                                                                                {
	#												print ALL ">gi\|10011$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
	#					        	        	                	print ALL "$var\n";
													print GTF "$coordis[0]\t10011$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$var\";\n";
													print GTF "$coordis[0]\t10011$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$var\";\n";
												}
												if ($strandis eq "-")
		                                                                                {
													$var_var = $var;
						        		        	                $rev_var = reverse $var_var;	
													if ($re_seq_size == 1)
													{
														if ($count == 1 || $count ==3)
														{
															$new_genomic_start=$new_genomic_start+1;$new_genomic_end=$new_genomic_end+1;
	#														print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
															print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
															print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
														}
														elsif ($count == 2)
														{
															if (($new_genomic_end-$jun_exons_cors[0]-2) >=0) 
															{
																$new_genomic_start=$new_genomic_start-2;$new_genomic_end=$new_genomic_end-2;
	#															print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
																print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
															}
															else
															{
																$new_genomic_reminder=$new_genomic_end-$jun_exons_cors[0]-2;
																if ($new_genomic_reminder == -2)
																{
																	$new_genomic_start=$new_genomic_start-2;$new_genomic_end=$exons_cors[1]-1;
																	print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																}
																else
																{
																	$new_genomic_start=$new_genomic_start-2;$new_genomic_end=$exons_cors[1];
																	print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																}
	#															print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
															}
														}
	#					        		        	                	print ALL "$rev_var\n";
													}
													elsif ($re_seq_size == 2)
													{
														if ($count ==1)
														{
															$new_genomic_start=$new_genomic_start+2;$new_genomic_end=$new_genomic_end+2;
	#														print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
															print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
															print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
														}
														elsif ($count == 2 || $count ==3)
														{
															if (($new_genomic_end-$jun_exons_cors[0]-1)>=0)
															{
																$new_genomic_start=$new_genomic_start-1;$new_genomic_end=$new_genomic_end-1;
	#															print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
																print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
															}
															else
															{
																$new_genomic_start=$new_genomic_start-1;$new_genomic_end=$exons_cors[1];
	#															print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
																print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
															}
														
														}
	#					        		        	                	print ALL "$rev_var\n";
													}
													else
													{
	#													print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
	#						        		        	                print ALL "$rev_var\n";
														print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
														print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
													}
												}
												$flag =1;last;
											}
										}
									}
								}
								else
								{
									$start = (($start+$exons_cors[0])-$exons_cors[1])-1;
								}if ($flag ==1) { last;}
							}
						}
					}
			}
	#}
						##################################################end of coordinate fetcher
					}
				}
				else
				{
					#print "\nstart tryptic Inside\n";
					#system ("date");		

					#coordinate_fetcher($var,$giaccession_new,$a,\%gene,$count,$seq_size,$for_pep_pos);
					####################################################start fo coordinate fetcher
					
#sub coordinate_fetcher
#	{
#			my ($var,$giaccession,$a,$hashref1,$count,$seq_size,$for_pep_pos)=@_;
			#print "\n<---- $count ---->\n";
			@pos=();
	#system ("date");
#			my %gene=%{$hashref1};
	#system ("date");	
			$seq_len = length ($var);
			#if ($seq_len >= 7)
			if ($seq_len >= 6 && $seq_len <= 35)
			{
					$giaccession++;
					$giaccession_1m = $giaccession;
					#print "$giaccession_1m\n";
					@coordis = split (/\|/, $gene{$a});
					$strandis = chop ($coordis[0]);
					if ($strandis eq "-")
					{
						$temp_pep_pos = reverse $for_pep_pos;
						$var = reverse $var;
						$result = index($temp_pep_pos, $var);
						$start = (($result+1)*3)-3; #((index is 0 based therefore +1)it is multiplied by 3 to get the triplet count and ) substracted by 2 because if triplet of M is ATG the start of the sequence will 3-2 i.e., ATG - TG = A;
						if ($count == 2)
						{
							$start = $start+2;
						}
						elsif ($count == 3)
						{
							$start = $start+1;
						}
					}
					else
					{
						$result = index($for_pep_pos, $var);
						$start = (($result+1)*3)-3; #((index is 0 based therefore +1)it is multiplied by 3 to get the triplet count and ) substracted by 2 because if triplet of M is ATG the start of the sequence will 3-2 i.e., ATG - TG = A;
						if ($count == 2)
						{
							$start = $start+1;
						}
						elsif ($count == 3)
						{
							$start = $start+2;
						}
					}
					$end = ($start+($seq_len*3))-1;
					$junc_end = $end - $start; $end_bckup = $end;
					@coordis = split (/\|/, $gene{$a});
					$strandis = chop ($coordis[0]);
					if ($coordis[1] !~ /\,/)
					{
						@exons_cors = split (/\-/, $coordis[1]);
						if ($strandis eq "+" || $strandis eq "\.")
						{
							$genomic_start = $exons_cors[0]+$start;
							$genomic_end = $genomic_start+($end-$start);
	#						print ALL ">gi\|10101$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$genomic_start\-$genomic_end\n";
        #	        		                print ALL "$var\n";
							print GTF "$coordis[0]\t10101$giaccession_1m\tCDS\t$genomic_start\t$genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$var\";\n";
						}
						elsif ($strandis eq "-")
						{
							$re_seq_size = $seq_size%3;
							if ($re_seq_size == 1)
							{
								if($count == 2)
								{		
									$genomic_start= ($exons_cors[0]+$start+1)-3;
									$genomic_end = $genomic_start+($end-$start);
								}
								else
								{
									$genomic_start= $exons_cors[0]+$start+1;
									$genomic_end = $genomic_start+($end-$start);
								}	
							}
							elsif ($re_seq_size == 2)
							{
								if($count == 2)
								{
									#print "santosh\n";
									$genomic_start= ($exons_cors[0]+$start+2)-3;
									$genomic_end = $genomic_start+($end-$start);
								}
								elsif($count == 3)
								{		
									$genomic_start= ($exons_cors[0]+$start+2)-3;
									$genomic_end = $genomic_start+($end-$start);
								}
								else
								{
									$genomic_start= $exons_cors[0]+$start+2;
									$genomic_end = $genomic_start+($end-$start);
								}	
							}
							else
							{
								$genomic_start= $exons_cors[0]+$start;
								$genomic_end = $genomic_start+($end-$start);
							}
							#print "behera\n";
	#						print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$genomic_start\-$genomic_end\n";
							$var_var = $var;
        		        	                $rev_var = reverse $var_var;
        #		        	                print ALL "$rev_var\n";
							print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$genomic_start\t$genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
						}
					}
					else 
					{
						@comma = split (/\,/, $coordis[1]);
						$start_bckup = $start;
						$index_comma=0;
						if ($strandis eq "+" || $strandis eq "-")
						{
							$flag=0;
							foreach $exonpair (@comma)
							{
								$index_comma++;
								@exons_cors = split (/\-/, $exonpair);
								if (($exons_cors[0]+$start) < $exons_cors[1])
								{
									$new_genomic_start = $exons_cors[0]+$start;
									$new_genomic_end = $new_genomic_start+($end-$start_bckup);
									if ($new_genomic_end <= $exons_cors[1])
									{
										if ($strandis eq "+")
										{
	#										print ALL ">gi\|10101$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
	#			        	        		                print ALL "$var\n";
											print GTF "$coordis[0]\t10101$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$var\";\n";
										}
										elsif ($strandis eq "-")
										{
											$re_seq_size = $seq_size%3;
											if ($re_seq_size == 1)
											{
												if($count == 2)
												{
													$new_genomic_start=($new_genomic_start+1)-3;
													$new_genomic_end=($new_genomic_end+1)-3;
												}
												else
												{
													$new_genomic_start=$new_genomic_start+1;
													$new_genomic_end=$new_genomic_end+1;
												}		
	#											print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
												$var_var = $var;
					        		        	                $rev_var = reverse $var_var;	
	#				        		        	                print ALL "$rev_var\n";
												print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
											}
											elsif ($re_seq_size == 2)
											{
												if($count == 2)
												{
													$new_genomic_start=($new_genomic_start+2)-3;
													$new_genomic_end=($new_genomic_end+2)-3;
												}
												elsif($count == 3)
												{
													$new_genomic_start=($new_genomic_start+2)-3;
													$new_genomic_end=($new_genomic_end+1)-2;
												}
												else
												{
													$new_genomic_start=$new_genomic_start+2;
													$new_genomic_end=$new_genomic_end+2;
												}	
												$var_var = $var;
					        		        	                $rev_var = reverse $var_var;	
	#											print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
	#				        		        	                print ALL "$rev_var\n";
												print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
											}
											else
											{
	#											print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
												$var_var = $var;
					        		        	                $rev_var = reverse $var_var;	
	#				        		        	                print ALL "$rev_var\n";
												print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
											}
										}
										last;
									}
									else
									{
										$p = $junc_end - ($exons_cors[1]-$new_genomic_start)-1;
										for ($i = $index_comma; $i <= $#comma; $i++)
										{
											@jun_exons_cors = split (/\-/, $comma[$i]);
											if (($jun_exons_cors[0]+$p) <= $jun_exons_cors[1])
											{
												#JUNCTIONS
												$new_genomic_end = $jun_exons_cors[0]+$p;
												if ($strandis eq "+")
		                                                                                {
	#												print ALL ">gi\|10011$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
	#					        	        	                	print ALL "$var\n";
													print GTF "$coordis[0]\t10011$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$var\";\n";
													print GTF "$coordis[0]\t10011$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$var\";\n";
												}
												if ($strandis eq "-")
		                                                                                {
													$var_var = $var;
						        		        	                $rev_var = reverse $var_var;	
													if ($re_seq_size == 1)
													{
														if ($count == 1 || $count ==3)
														{
															$new_genomic_start=$new_genomic_start+1;$new_genomic_end=$new_genomic_end+1;
	#														print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
															print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
															print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
														}
														elsif ($count == 2)
														{
															if (($new_genomic_end-$jun_exons_cors[0]-2) >=0) 
															{
																$new_genomic_start=$new_genomic_start-2;$new_genomic_end=$new_genomic_end-2;
	#															print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
																print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
															}
															else
															{
																$new_genomic_reminder=$new_genomic_end-$jun_exons_cors[0]-2;
																if ($new_genomic_reminder == -2)
																{
																	$new_genomic_start=$new_genomic_start-2;$new_genomic_end=$exons_cors[1]-1;
																	print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																}
																else
																{
																	$new_genomic_start=$new_genomic_start-2;$new_genomic_end=$exons_cors[1];
																	print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																}
	#															print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
															}
														}
	#					        		        	                	print ALL "$rev_var\n";
													}
													elsif ($re_seq_size == 2)
													{
														if ($count ==1)
														{
															$new_genomic_start=$new_genomic_start+2;$new_genomic_end=$new_genomic_end+2;
	#														print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
															print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
															print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
														}
														elsif ($count == 2 || $count ==3)
														{
															if (($new_genomic_end-$jun_exons_cors[0]-1)>=0)
															{
																$new_genomic_start=$new_genomic_start-1;$new_genomic_end=$new_genomic_end-1;
	#															print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
																print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
															}
															else
															{
																$new_genomic_start=$new_genomic_start-1;$new_genomic_end=$exons_cors[1];
	#															print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
																print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
															}
														
														}
	#					        		        	                	print ALL "$rev_var\n";
													}
													else
													{
	#													print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
	#						        		        	                print ALL "$rev_var\n";
														print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
														print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
													}
												}
												$flag =1;last;
											}
										}
									}
								}
								else
								{
									$start = (($start+$exons_cors[0])-$exons_cors[1])-1;
								}if ($flag ==1) { last;}
							}
						}
					}
			}
	#}
					####################################################end of cooriante fetcher
					#print "\nEnd tryptic Inside\n";
					#system ("date");		

				}
			}
		}	
	#}
		###############
	#print "\nend tryptic\n";
	#system ("date");		
		
		$protein="";
		@split_seq="";

	}
#}
	$mycount++;
	#if ($mycount%1000==0){print "\t$mycount\n";
	#delete $seq{$a};

	#system ("date");}
	
}
#Create FASTA file for sequence between stop to stop. 
stop2stop(\%codon,\%seq_s2s,\%more);

############################################################################################################################

#####
#subroutine
#######
################
##
#############

sub unique{
	my ($arrref)=@_;
	my @array=@{$arrref};
	my %uniqueHash="";
	foreach (@array)
	{
		$uniqueHash{$_}.="";
	}
	my @unique=keys %uniqueHash;
	return \@unique;
}
	
###################
#Enzymes#####
#############

sub trypsin{

	my ($seq)=@_;                
	$seq =~ s/([RK])/$1=/g;
	$seq =~ s/=P/P/g;
	my @split_seq = split (/=/, $seq);
	return \@split_seq;
}
sub argC{

	my ($seq)=@_;              
	$seq =~ s/(R)/$1=/g;
	#$seq =~ s/=P/P/g;
	my @split_seq = split (/=/, $seq);
	return \@split_seq;
}
sub lysC{
	my ($seq)=@_;
	$seq =~ s/(L)/$1=/g;
	my @split_seq = split (/=/,$seq);
	return \@split_seq;
}
sub lysN{
	my ($seq)=@_;
	$seq =~ s/(L)/=$1/g;
	my @split_seq = split (/=/,$seq);
	return \@split_seq;
}
sub aspN{
	my ($seq)=@_;
	$seq =~ s/(D)/=$1/g;
	my @split_seq = split (/=/,$seq);
	return \@split_seq;
}
sub gluC{
	my ($seq)=@_;
	$seq =~ s/([ED])/$1=/g;
	$seq =~ s/=P/P/g;
	my @split_seq = split (/=/,$seq);
	return \@split_seq;
}
sub chymotrypsin{

	my ($seq)=@_;  
	$seq =~ s/WM/XX/g;              
	$seq =~ s/([WYF])/$1=/g;
	$seq =~ s/XX/WM/g;
	$seq =~ s/=P/P/g;
	my @split_seq = split (/=/, $seq);
	return \@split_seq;
}


#########################
## MC 1
###############
sub mc1{
	my ($arrref1)=@_;
	my @split_seq=@{$arrref1};
	my $split_seq_size=scalar(@split_seq);
	my @seq_split_mc="";
	for ($mc =0; $mc <= $split_seq_size; $mc++)
	{
		$miscleavage = "$split_seq[$mc]"."$split_seq[$mc+1]";
		unless($miscleavage=~m/\*/){push @seq_split_mc, $miscleavage;}
	}
	
	return \@seq_split_mc;
}
#########################
## MC 2
#################
sub mc2{
	my ($arrref1)=@_;
	my @split_seq=@{$arrref1};
	my $split_seq_size=scalar(@split_seq);
	my @seq_split_mc="";
	for ($mc =0; $mc <= $split_seq_size; $mc++)
	{
		$miscleavage = "$split_seq[$mc]"."$split_seq[$mc+1]"."$split_seq[$mc+2]";
		unless($miscleavage=~m/\*/){push @seq_split_mc, $miscleavage;}
	}
	
	return \@seq_split_mc;
}



#########
#subroutine2
###########

###########
#coordinate-fetcher
############

#############
#subroutine 3####
#############


sub stop2stop{

($hashref1,$hashref2,$hashref3)=@_;

my %codon= %{$hashref1};
my %seq_s2s= %{$hashref2};
my %more= %{$hashref3};




$giaccession=0;
$intnum = 1; 



while (($a,$b) = each %seq_s2s) #hash with transcript id as key and nucleotide sequence as value
{
	
	@sequence = split (//, $b); #$a is trancript id and $b is nucleotide
	$seq_size = @sequence;
	for ($j =0; $j<=2; $j++)
	{
		for ($i=$j; $i<= $seq_size; $i=$i+3)
		{
			$triplet = "$sequence[$i]"."$sequence[$i+1]"."$sequence[$i+2]";#taking codons
			$protein .= "$codon{$triplet}"; #this hash-codon should be taken as input when calling this subroutine
		}
		$prt_len = length $protein;#length of protein
		$seq = $protein; 
		$count = $j+1;
		@stars = split (/\*/, $protein);#split with stop codons
		$sub_count=0;
		foreach $xs (@stars)
		{
			$xs_len = length ($xs);
			$sub_count++;
			if ($xs_len >= 7)
			{
				$new_temp  =  $xs;
				$new_temp =~ s/(K|R)/$1=/g;
				$new_temp =~ s/=P/P/g;
				$new_temp =~ s/(K|R)=/#/g;
				$temp_v = 0; $char="#";$flag=0;
				my $offset = 0;
				my $result = index($new_temp, $char, $offset);#return the position of $char in $new_temp after at/after $offset i.e 0
				if ($result != -1) # result = -1 when the $char didn't match $new_temp
				{
					while ($result != -1) #search until it matches #
					{
						$new_v = $temp_v;
						$temp_v = $result+1;
						$diff = $temp_v-$new_v;
						if ($diff >= 7) {$flag =1; }
						$offset = $result + 1; #incrementing the offset i.e the starting of search
						$result = index($new_temp, $char, $offset);
					}
					if ($flag ==1)
					{
						$intlen = length ($intnum);
						if ($intlen == 1) {$someGiNum = "1000000$intnum";}
						if ($intlen == 2) {$someGiNum = "100000$intnum";}
						if ($intlen == 3) {$someGiNum = "10000$intnum";}
						if ($intlen == 4) {$someGiNum = "1000$intnum";}
						if ($intlen == 5) {$someGiNum = "100$intnum";}
						if ($intlen == 6) {$someGiNum = "10$intnum";}
						if ($intlen == 7) {$someGiNum = "1$intnum";}
						
						print OUT ">gi|$someGiNum|ref|$a| $a\#$more{$a}\#PF$count\_PC$sub_count\n$xs\n";#hash more to be taken
						$intnum++;
						$result=0;$diff=0;$new_v=0;
					}
				}
			}
		}
		$protein="";
	}
}

}
