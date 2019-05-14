#!/usr/bin/perl
################################
### WIthout Subroutine  Extension to new5.pl###
################################



unless(scalar (@ARGV) >= 6 )
{
	print "Give the proper format: perl GenoART.pl 0/1/2_MC Enzyme_name input.fasta genetic_code proteotypic_peptide_Length ProteinLength\n";
	exit;
}
else
{
	$infile=$ARGV[2];
	$cleavage=$ARGV[0];
	$enzyme=$ARGV[1];
	$genetic_code_tt = $ARGV[3];
	$minProPepLength = $ARGV[4];
	$minProteinLength = $ARGV[5];
	$outputSuffixName = $ARGV[6];
}

$out_m2s=$out_s2s=$out_gtf=$in_file=$infile;
($type)=(split /\./,$infile)[-1];
$in_file=~s/\.$type/\_op\.$type/;

$out_gtf=~s/\.$type/\_pep\_custom\_$outputSuffixName\.gtf/;
$out_s2s=~s/\.$type/\_pep\_custom\_$outputSuffixName\.$type/;
$out_m2s=~s/\.$type/\_pep\_custom\_m2s\_$outputSuffixName\.$type/;

open FILE, "$in_file" or die $!;
open (GENETIC, "bin/genetic_code/$genetic_code_tt") or die $!;

while ($genetic = <GENETIC>)
{
	chomp($genetic);
	$genetic =~ s/\r//g;
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
#open (ALL, ">$out_fasta") or die $!;
open (GTF, ">$out_gtf") or die $!;
open (OUT, ">$out_s2s") or die $!;
open (M2SOUT, ">$out_m2s") or die $!;
#################################

#Create FASTA file for sequence between stop to stop. 
stop2stop(\%codon,\%seq_s2s,\%more);

print "\nCreating 6/3 - frame database is complete\n\n";

print "Creating custom GTF is in progress\n\n";

while (($a,$b) = each %seq)  
{
	@sequence = split (//, $b);
	$seq_size = @sequence; 
	my @from_mc1=@seq_split_mc=();
	$seq_size=scalar(@sequence);
	for ($j =0; $j<=2; $j++)       # Creatign three frames j= 0,1,2
	{
		for ($i=$j; $i<= $seq_size; $i=$i+3) # take the three ltters and compare it with codon
		{
			$triplet = "$sequence[$i]"."$sequence[$i+1]"."$sequence[$i+2]"; #mergin the three letter for one codon
			if (exists $codon{$triplet})
			{
				$protein .= "$codon{$triplet}"; # Fetch the amino acids for the above three letter if exists 
			}
			else
			{
					$protein .= "*"; # Insert star in place of amino acids where the above three letters are NNN or NNA or NNT or NNG or NNC or other combinations of xNx
			}
		}
		$prt_len = length $protein; # taking the length of concatinated amino acids
		$seq = $protein; #Storing the ptoeins in seq
		#print "------**>$seq\n";
		$for_pep_pos = $protein;
		#print "$a\n$seq\n";
		$count = $j+1; 
		if($enzyme eq "trypsin")
		{
			$arrref=trypsin($seq);
			@split_seq=@{$arrref};
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
			$ref_from_mc1=mc1(\@split_seq);
			@from_mc1=@{$ref_from_mc1};
			push @from_mc1,@split_seq;
			$arrref1=mc2(\@split_seq);
			@seq_split_mc=@{$arrref1};
			push @seq_split_mc,@from_mc1;
			@split_seq=@seq_split_mc;@seq_split_mc=@from_mc1="";

		}
		$uref=unique(\@split_seq);
		@split_seq=@{$uref};
		#print "-->@split_seq\n\n";
		##############
		foreach $var1 (@split_seq)
		{
			@spStars = ();
			if ($var1 =~ /\*/) # if the tryptic peptide contains *
			{	
				@spStars = split (/\*/, $var1);
				
			}
			else
			{
				@spStars = $var1;
			}
			foreach $var (@spStars)
			{
			
				@mysearchArraySeq = ();@mets = ();
				$giaccession_new++;
				####################################################start fo coordinate fetcher
				@pos=();
				$seq_len = length ($var);
				if ($seq_len >= $minProPepLength)
				#if ($seq_len >= 6 && $seq_len <= 35)
				{
					#print "-->$var\n";
					$peptide_seq_M = $var;
					@mysearchArraySeq = ("$peptide_seq_M");
					#print "--> $peptide_seq_M\n";
					if ($peptide_seq_M =~ /M/)
					{
				         	$newmetSeq = $peptide_seq_M;
				                $newmetSeq =~ s/M/\#M/;
					        @mets = split (/\#/, $newmetSeq);
					        $theMString = $mets[-1];
					        foreach $sn (0..2)
				        	{
			        	        	push @mysearchArraySeq, substr($theMString, $sn,);
					        }
					}
					#print "-------> @mysearchArraySeq\n";
					foreach $peptide_seq (@mysearchArraySeq)
					{
						#if ($a =~ /ENST00000480537.5/) {print "I am in this transcript now: ENST00000480537.5\n";}
						#if ($peptide_seq =~ /PVAMASYYEILDVPR/) { print " I am here: $peptide_seq\n";} 
						$result=$start=$end=$junc_end=$genomic_start=$genomic_end=0; 
						$seq_len = length($peptide_seq);
						if ($seq_len >= $minProPepLength)
						{
							#print "$peptide_seq\n";
							$var = $peptide_seq;
							$var_str =$var;
							$giaccession++;
							$giaccession_1m = $giaccession;
							#print "$giaccession_1m\n";
							@coordis = split (/\|/, $gene{$a});
							$strandis = chop ($coordis[0]);
							if ($strandis eq "-")
							{
								$temp_pep_pos = reverse $for_pep_pos;
								$var = reverse $var;
								#print "$temp_pep_pos\n$peptide_seq\t-->$var\n";
								$result = index($temp_pep_pos, $var);
								#print "$peptide_seq\t --> $result\n";
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
		        #		        	print ALL "$rev_var\n";
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
				#													print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
																	print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																	print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																}
																elsif ($count == 2)
																{
																	if (($new_genomic_end-$jun_exons_cors[0]-2) >=0) 
																	{
																		$new_genomic_start=$new_genomic_start-2;$new_genomic_end=$new_genomic_end-2;
				#														print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";
																		print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																		print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";															}
																	else
																	{
																		$new_genomic_reminder=$new_genomic_end-$jun_exons_cors[0]-2;
																		if ($new_genomic_reminder == -2)
																		{
																			$new_genomic_start=$new_genomic_start-2;$new_genomic_end=$exons_cors[1]-1;
																			print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";														}
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
			#															print ALL ">gi\|10010$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$exons_cors[1],$jun_exons_cors[0]\-$new_genomic_end\n";																	print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$new_genomic_start\t$exons_cors[1]\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																		print GTF "$coordis[0]\t10010$giaccession_1m\tCDS\t$jun_exons_cors[0]\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";
																	}
																	else
																	{
																		$new_genomic_start=$new_genomic_start-1;$new_genomic_end=$exons_cors[1];
			#															print ALL ">gi\|10110$giaccession_1m\|ref|$a\| $a\#$strandis$count\#$coordis[0]\:$new_genomic_start\-$new_genomic_end\n";
																		print GTF "$coordis[0]\t10110$giaccession_1m\tCDS\t$new_genomic_start\t$new_genomic_end\t\.\t$strandis\t\.\tgene_id \"$a\"\; transcript_id \"$rev_var\";\n";															}
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
					}
				}
			####################################################end of cooriante fetcher
			}
		}	
		$protein="";
		@split_seq="";
	}
	$mycount++;
}
#Create FASTA file for sequence between stop to stop. 
#stop2stop(\%codon,\%seq_s2s,\%more);

############################################################################################################################

#####
#subroutine
#######
sub unique{
	my ($arrref)=@_;
	my @array=@{$arrref};
	my %uniqueHash=();
	foreach (@array)
	{
		$uniqueHash{$_}.="";
	}
	my @unique=keys %uniqueHash;
	return \@unique;
}
	
#############
#Enzymes#####
#############

sub trypsin{

	my ($seq)=@_;                
	$seq =~ s/([RK])/$1=/g;
	$seq =~ s/=P/P/g;
	my @split_seq = split (/=/, $seq);
	#print "$seq\n--> @split_seq\n";
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
	$seq =~ s/(K)/$1=/g;
	my @split_seq = split (/=/,$seq);
	return \@split_seq;
}
sub lysN{
	my ($seq)=@_;
	$seq =~ s/(K)/=$1/g;
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


###############
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
		#print "$miscleavage\n";
		#unless($miscleavage=~m/\*/){push @seq_split_mc, $miscleavage;}
		push @seq_split_mc, $miscleavage;
	}
	
	return \@seq_split_mc;
}
#################
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
		#unless($miscleavage=~m/\*/){push @seq_split_mc, $miscleavage;}
		push @seq_split_mc, $miscleavage;
	}
	
	return \@seq_split_mc;
}

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
				if ($xs_len >= $minProteinLength)
				{
					$intlen = length ($intnum);
					$someGiNum = "10$intnum";
					if ($xs =~ /M/)
                   	 		{
						if ($xs =~ /^M/)
						{
							if (!exists $checkRedundancyMet{$xs})
							{
								$checkRedundancyMet{$xs}="";
								print OUT ">gi|$someGiNum|ref|$a| \n$xs\n";#hash more to be taken #Date modified: 21/10/2018
								print M2SOUT ">gi|$someGiNum|ref|$a| \n$xs\n";
								$intnum++;
							}
						}
						else
						{
							$NewMetSeq = $xs;
							$NewMetSeq =~ s/M/\#M/;
							@NewMets = split (/\#/, $NewMetSeq);
							$theNewMString = $NewMets[-1];
							if ($theNewMString >=$minProteinLength)
							{
								if (!exists $checkRedundancyMet{$xs})
								{
									$checkRedundancyMet{$xs}="";
									print OUT ">gi|$someGiNum|ref|$a| \n$xs\n";#hash more to be taken #Date modified: 21/10/2018
									print M2SOUT ">gi|$someGiNum|ref|$a| \n$xs\n";
									$intnum++;
								}
							}
						}
                    			}
					else
					{
						if (!exists $checkRedundancy{$xs})
						{
							$checkRedundancy{$xs}="";
							print OUT ">gi|$someGiNum|ref|$a| \n$xs\n";#hash more to be taken #Date modified: 21/10/2018
							$intnum++;
						}
					}
					#print OUT ">gi|$someGiNum|ref|$a| \n$xs\n";#hash more to be taken #Date modified: 21/10/2018
					#print OUT ">gi|$someGiNum|ref|$a| $a\#$more{$a}\#PF$count\_PC$sub_count\n$xs\n";#hash more to be taken
					$result=0;$diff=0;$new_v=0;
				}
			}
			$protein="";
		}
	}
}
