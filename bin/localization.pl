#!/usr/bin/perl
$gtf_input = "$ARGV[0]"; #Input GTF for categorization
$fasta_input = "$ARGV[1]"; #Input protein FASTA file 
$ref_annotation = "$ARGV[2]"; #Input reference annotation GTF for comparision
$input_gssps =  "$ARGV[3]"; #Input GSSPs file after proteogenomics analysis
@dots_out = split (/\./, $input_gssps);
$output = "$dots_out[0]"."_categorized.txt";
$output_gtf = "$dots_out[0]"."_categorized.gtf";
$output_prt = "$dots_out[0]"."_prototeinCoding.txt";
$notmapped_output = "$dots_out[0]"."_notMapped.txt";
$output_gtf_prt = "$dots_out[0]"."_proteinCoding.gtf";
open (INP, "$gtf_input") or  die $!;
open (PRT, "$fasta_input") or die $!;
while ($prt_seq = <PRT>)
{
	chomp ($prt_seq);
	$prt_seq =~ s/\r//g;
	if ($prt_seq !~ /^>/)
	{
		$protein_coding_seq .= "$prt_seq";
	}	
}
print "Reading the protein sequence is complete\n";
open (GTF, "$ref_annotation") or die $!;
while (<GTF>)
{
	chomp ($_);
	if ($_ !~ /^#/)
	{
		@cols = split (/\t/, $_);
		if ($cols[2] !~ /[G|g]ene/ && $cols[2] !~ /[E|e]xon/ && $cols[2] !~ /[C|c][D|d][S|s]/ && $cols[2] !~ /[U|u][T|t][R|r]/ && $cols[2] !~ /[S|s]tart_codon/ && $cols[2] !~ /[S|s]tart_codon/ && $cols[2] !~ /[S|s]top_codon/ && $cols[2] !~ /[M|m]atch/ && $cols[2] !~ /[C|c]entromere/ && $cols[2] !~ /[C|c]ontig/ && $cols[2] !~ /[D|d]_[L|l]oop/ && $cols[2] !~ /[E|e]nhancer/ && $cols[2] !~ /[L|l]ow_[C|c]omplexity_[R|r]egion/ && $cols[2] !~ /[P|p]romoter/ && $cols[2] !~ /[R|r]egion/ && $cols[2] !~ /[S|s]elenocysteine/ && $cols[2] !~ /[S|s]equence_[F|f]eature/ && $cols[2] !~ /[T|t]andem_[R|r]epeat/)
		{
			$main_key = "$cols[0]\#$cols[6]"; $trans_cords = "$cols[3]\#$cols[4]";
			$transcript{$main_key} .= "$trans_cords;";
			@desc_cols = split (/\;/, $cols[8]);
			$gene_type_c= "-";$gene_name_c ="-";$description_c = "-";
			foreach $dtype (@desc_cols)
			{
				if ($dtype =~ /[G|g]ene[\:|\=]/ || $dtype =~ /[P|p]arent\=/ || $dtype =~ /[G|g]ene_[N|n]ame/)
				{
					$gene_name_c = $dtype; $gene_name_c =~ s/.*[G|g]ene[\:|\=]//;$gene_name_c =~ s/.*[P|p]arent\=//; $gene_name_c =~ s/.*[G|g]ene_[N|n]ame//;$gene_name_c =~ s/ //g; $gene_name_c =~ s/"//g;
				}
				if ($dtype =~ /[P|p]roduct\=/ || $dtype =~ /[D|d]escription[\=|\:]/)
				{
					$description_c = $dtype; $description_c =~ s/.*[P|p]roduct\=//; $description_c =~ s/.*[D|d]escription[\=|\:]//; $description_c =~ s/"//g; 
				}
				if ($dtype =~ /[G|g]ene_type/)
				{
					$gene_type_c = "$dtype"; $gene_type_c=~ s/.*[G|g]ene_type//;$gene_type_c=~ s/"//g;
				}
			}
			$transcript_type{$trans_cords} = "$gene_name_c\t$description_c\t$gene_type_c";
		}
		if ($cols[2] =~ /[E|e]xon/)
		{
			$ecords = "Exon:$trans_cords"; $exon_cords = "$cols[3]\#$cols[4]";
			$exon_key{$ecords} .= "$exon_cords;";
		}
		if ($cols[2] =~ /[C|c][D|d][S|s]/)
		{
			$ccords = "CDS:$trans_cords"; $CDS_cords = "$cols[3]\#$cols[4]";
			$cds_key{$ccords} .= "$CDS_cords;";
		}
		if ($cols[2] =~ /[S|s]tart_codon/)
		{
			$startcords = "START:$trans_cords"; $start_cords = "$cols[3]\#$cols[4]";
			$start_key{$startcords} .= "$start_cords;";
		}
		if ($cols[2] =~ /[S|s]top_codon/)
		{
			$endcords = "STOP:$trans_cords"; $end_cords = "$cols[3]\#$cols[4]";
			$end_key{$endcords} .= "$end_cords;";
		}
	}
}
print "Reading GTF is complete!!!\n";
open (OUT, ">$output") or die $!;
open (GTF_OUT, ">$output_gtf") or die $!;
open (PRT_OUT, ">$output_prt") or die $!;
open (PRT_GTF_OUT, ">$output_gtf_prt") or die $!;
open (INP_GSSP, "$input_gssps") or die $!;
print "$input_gssps\n";
while ($ingssp = <INP_GSSP>)
{
	chomp ($ingssp);
	$ingssp =~ s/\r//g;
	$ingssp =~ s/\"//g;
	@gsspeps = split (/\t/, $ingssp);
	if ($ingssp =~ /[S|s]equence/)
	{
		for ($hi =0; $hi <= $#gsspeps; $hi++)
		{
			if ($gsspeps[$hi] =~ /^[P|p]eptide/ || $gsspeps[$hi] =~ /^[S|s]equence/)
			{
				$peptide_element = $hi;
				last;
			}
		}
	}
	else
	{
		if ($gsspeps[$peptide_element] =~ /\./)
		{
			@dots_gssps = split (/\./, $gsspeps[$peptide_element]);
			$dots_gssps[1] =~ tr/[a-z]/[A-Z]/;
			$req_peptide_seq{$dots_gssps[1]} .= "$ingssp*";
		}
		else
		{
			$gsspeps[$peptide_element] =~ tr/[a-z]/[A-Z]/;
			$req_peptide_seq{$gsspeps[$peptide_element]} .= "$ingssp*";
		}
	}
}

@gssp_size = keys %req_peptide_seq;
$size_gssp = @gssp_size;

print "Reading input GSSP is complete!!!\n";
print "\nLocalizing and categorizing the peptides!!!\n";
$count=0;
print "Total number of peptides in the GSSP list is: $size_gssp\n";
while ($geno = <INP>)
{
	@mysearchArraySeq = ();
	$novel_transcript = 1;
	chomp ($geno);
	$geno =~ s/\r//g;
	@gsplit = split (/\t/, $geno);
	$strand = "$gsplit[6]";
	$specific_chr = "$gsplit[0]\#$gsplit[6]";
	@pep_coln = split (/\"/, $gsplit[8]);
	$peptide_seq = "$pep_coln[3]";
	if (exists $req_peptide_seq{$peptide_seq})
	{
		$mapped_peptides{$peptide_seq} = "mapped";
		print "Reading peptide\: $peptide_seq \n";
		if ($gsplit[1] =~ /^10101/) {$splice = "Exonic";}
		if ($gsplit[1] =~ /^10011/) {$splice = "Exon junction";}
		if ($gsplit[1] =~ /^10110/) {$splice = "Exonic";}
		if ($gsplit[1] =~ /^10010/) {$splice = "Exon junction";}
		if (exists $transcript{$specific_chr}) #Chromosome and Strand should match;
		{
			@t_coordinates = split (/\;/, $transcript{$specific_chr});
			@cats = "";
			foreach $tc (@t_coordinates)
			{
				$stcodon = "START:$tc";$encodon = "STOP:$tc";
				$cdsst = "CDS:$tc"; $utrst = "UTR:$tc"; $exonst = "Exon:$tc";
				@hash_co = split (/\#/, $tc);
				if ($strand eq "\+")
				{
					if ($gsplit[3] >= $hash_co[0] && $gsplit[4] <= $hash_co[1] )#With in transcript;
					{
						$novel_transcript = 0;
						$TYPE = $transcript_type{$tc};
						if (exists $start_key{$stcodon}) #If exists start codons;
						{
							@iniators = split (/\;/, $start_key{$stcodon});
							@start_cdn = split (/\#/, $iniators[0]);
							if ($gsplit[3] < $start_cdn[0]) {$category = "5' UTR";$novel_exon=0; push @cats, $category;} 
							if ($gsplit[3] < $start_cdn[0] && $gsplit[4] >= $start_cdn[1]) 
							{
								$category = "N_terminal protein extension";
								$novel_exon=0;
								push @cats, $category;
							}
						}
						if (exists $end_key{$encodon}) #If exists end/stop codons;
						{
							@terminators = split (/\;/, $end_key{$encodon});
							@end_cdn = split (/\#/, $terminators[0]);
							if ($gsplit[4] > $end_cdn[1]) {$category = "3' UTR";$novel_exon=0; push @cats, $category; } 
							if ($gsplit[3] <= $end_cdn[0] && $gsplit[4] > $end_cdn[1]) {$category = "C_terminal protein extension"; $novel_exon=0;push @cats, $category; } 
						}
						if (exists $cds_key{$cdsst}) #If exists, CDS coordinates;
						{
							@cds_semi = split (/\;/, $cds_key{$cdsst});
							@cds_st_coors = split (/\#/, $cds_semi[0]); @cds_en_coors = split (/\#/, $cds_semi[-1]);
							if ($gsplit[4] < $cds_st_coors[0]) {$category = "5' UTR"; $novel_exon=0;push @cats, $category;}  
							if ($gsplit[3] > $cds_en_coors[1]) {$category = "3' UTR"; $novel_exon=0; push @cats, $category;}  
							if ($gsplit[3] < $cds_st_coors[0] && $gsplit[4] > $cds_st_coors[1]) {$category = "N_terminal protein extension";push @cats, $category; $novel_exon=0;}
							if ($gsplit[4] > $cds_en_coors[1] && $gsplit[3] < $cds_en_coors[0]) {$category = "C_terminal protein extension";push @cats, $category; $novel_exon=0;}
							if ($gsplit[3] >= $cds_st_coors[0] && $gsplit[4] <= $cds_en_coors[1]) #With in CDS;
							{
								$novel_exon =1;
								foreach $cs (@cds_semi)
								{
									@each_cds_coors = split (/\#/, $cs);
									if ($gsplit[3] >= $each_cds_coors[0] && $gsplit[4] <= $each_cds_coors[1])
									{
										$category = "Alternate frame of translation";
										$novel_exon=0; #last;
										push @cats, $category;
									}
									elsif ($gsplit[3] < $each_cds_coors[0] && $gsplit[4] >= $each_cds_coors[0] && $gsplit[4] <= $each_cds_coors[1])
									{
										$category = "5' Exon extension";
										$novel_exon=0; #last;
										push @cats, $category;
									}
									elsif ($gsplit[3] >= $each_cds_coors[0] && $gsplit[3] <= $each_cds_coors[1] && $gsplit[4] > $each_cds_coors[1])
									{
										$category = "3' Exon extension";
										$novel_exon=0; #last;
										push @cats, $category;
									}
								}
								if ($novel_exon == 1) {$category = "Novel Exon";push @cats, $category;} #last;
							}
						}
						if ($exon_key{$exonst}) #If exists, in Exon
						{
							@exon_semi = split (/\;/, $exon_key{$exonst});
							$novel_exon =1;
							foreach $es (@exon_semi)
							{
								@each_exon_coors = split (/\#/, $es);
								if ($gsplit[3] >= $each_exon_coors[0] && $gsplit[4] <= $each_exon_coors[1])
								{
										$category = "Alternate frame of translation";
										$novel_exon=0; #last;
										push @cats, $category;
								}
								elsif ($gsplit[3] < $each_exon_coors[0] && $gsplit[4] >= $each_exon_coors[0] && $gsplit[4] <= $each_exon_coors[1])
								{
										$category = "5' Exon extension";
										$novel_exon=0; #last;
										push @cats, $category;
								}
								elsif ($gsplit[3] >= $each_exon_coors[0] && $gsplit[3] <= $each_exon_coors[1] && $gsplit[4] > $each_exon_coors[1])
								{
										$category = "3' Exon extension";
										$novel_exon=0; #last;
										push @cats, $category;
								}
							}
							if ($novel_exon == 1) {$category = "Novel Exon";push @cats, $category;} 
						}
						$novel_transcript = 0;
					}
					elsif ($gsplit[3] < $hash_co[0] && $gsplit[4] > $hash_co[0] && $gsplit[4] <= $hash_co[1])
					{
						$category = "5' Gene extension";
						$novel_transcript = 0;
						push @cats, $category;
					}
					elsif ($gsplit[3] < $hash_co[1] && $gsplit[3] >= $hash_co[0] && $gsplit[4] > $hash_co[1])
					{
						$category = "3' Gene extension";
						$novel_transcript = 0;
						push @cats, $category;
					}
				}
				elsif ($strand eq "\-")
				{
					if ($gsplit[3] >= $hash_co[0] && $gsplit[4] <= $hash_co[1] )#With in transcript;
					{
						$novel_transcript = 0;
						$TYPE = $transcript_type{$tc};
						if (exists $start_key{$stcodon}) #If exists start codons;
						{
							@iniators = split (/\;/, $start_key{$stcodon});
							@start_cdn = split (/\#/, $iniators[0]);
							if ($gsplit[3] > $start_cdn[1]) {$category = "5' UTR";$novel_exon=0; push @cats, $category;} 
							if ($gsplit[3] <= $start_cdn[0] && $gsplit[4] > $start_cdn[1]) 
							{
								$category = "N_terminal protein extension";
								$novel_exon=0;
								push @cats, $category;
							}
						}
						if (exists $end_key{$encodon}) #If exists end/stop codons;
						{
							@terminators = split (/\;/, $end_key{$encodon});
							@end_cdn = split (/\#/, $terminators[0]);
							if ($gsplit[4] < $end_cdn[0]) {$category = "3' UTR";$novel_exon=0; push @cats, $category; } 
							if ($gsplit[3] < $end_cdn[0] && $gsplit[4] >= $end_cdn[1]) {$category = "C_terminal protein extension"; $novel_exon=0;push @cats, $category; } 
						}
						if (exists $cds_key{$cdsst}) #If exists, CDS coordinates;
						{
							@cds_semi = split (/\;/, $cds_key{$cdsst});
							@cds_st_coors = split (/\#/, $cds_semi[0]); @cds_en_coors = split (/\#/, $cds_semi[-1]);
							if ($cds_st_coors[0] >= $cds_en_coors[0]) 
							{
								@rev_cds = reverse @cds_semi;
								@cds_semi = @rev_cds;
								@cds_st_coors = split (/\#/, $cds_semi[0]); @cds_en_coors = split (/\#/, $cds_semi[-1]);
							}
							if ($gsplit[4] < $cds_st_coors[0]) {$category = "3' UTR"; $novel_exon=0; push @cats, $category;} # last;
							if ($gsplit[3] > $cds_en_coors[1]) {$category = "5' UTR";$novel_exon=0; push @cats, $category;} # last;
							if ($gsplit[3] < $cds_st_coors[0] && $gsplit[4] >= $cds_st_coors[0]) {$category = "C_terminal protein extension";$novel_exon=0; push @cats, $category;} # last;
							if ($gsplit[4] > $cds_en_coors[1] && $gsplit[3] <= $cds_en_coors[1]) {$category = "N_terminal protein extension";$novel_exon=0; push @cats, $category;} # last;
							if ($gsplit[3] >= $cds_st_coors[0] && $gsplit[4] <= $cds_en_coors[1]) #With in CDS;
							{
								$novel_exon =1;
								foreach $cs (@cds_semi)
								{
									@each_cds_coors = split (/\#/, $cs);
									if ($gsplit[3] >= $each_cds_coors[0] && $gsplit[4] <= $each_cds_coors[1])
									{
										$category = "Alternate frame of translation";
										$novel_exon=0; #last;
										push @cats, $category;
									}
									elsif ($gsplit[3] < $each_cds_coors[0] && $gsplit[4] >= $each_cds_coors[0] && $gsplit[4] <= $each_cds_coors[1])
									{
										$category = "3' Exon extension";
										$novel_exon=0; #last;
										push @cats, $category;
									}
									elsif ($gsplit[3] >= $each_cds_coors[0] && $gsplit[3] <= $each_cds_coors[1] && $gsplit[4] > $each_cds_coors[1])
									{
										$category = "5' Exon extension";
										$novel_exon=0; #last;
										push @cats, $category;
									}
								}
								if ($novel_exon == 1) {$category = "Novel Exon";push @cats, $category;} #last;
							}
						}
						if ($exon_key{$exonst}) #If exists, in Exon
						{
							@exon_semi = split (/\;/, $exon_key{$exonst});
							$novel_exon =1;
							foreach $es (@exon_semi)
							{
								@each_exon_coors = split (/\#/, $es);
								if ($gsplit[3] >= $each_exon_coors[0] && $gsplit[4] <= $each_exon_coors[1])
								{
										$category = "Alternate frame of translation";
										$novel_exon=0; #last;
										push @cats, $category;
								}
								elsif ($gsplit[3] < $each_exon_coors[0] && $gsplit[4] >= $each_exon_coors[0] && $gsplit[4] <= $each_exon_coors[1])
								{
										$category = "3' Exon extension";
										$novel_exon=0; #last;
										push @cats, $category;
								}
								elsif ($gsplit[3] >= $each_exon_coors[0] && $gsplit[3] <= $each_exon_coors[1] && $gsplit[4] > $each_exon_coors[1])
								{
										$category = "5' Exon extension";
										$novel_exon=0; #last;
										push @cats, $category;
								}
							}
							if ($novel_exon == 1) {$category = "Novel Exon";push @cats, $category;} #last;
						}
						$novel_transcript = 0;
					}
					elsif ($gsplit[3] < $hash_co[0] && $gsplit[4] > $hash_co[0] && $gsplit[4] <= $hash_co[1])
					{
						$category = "3' Gene extension";
						$novel_transcript = 0; #last;
						push @cats, $category;
					}
					elsif ($gsplit[3] < $hash_co[1] && $gsplit[3] >= $hash_co[0] && $gsplit[4] > $hash_co[1])
					{
						$category = "5' Gene extension";
						$novel_transcript = 0; #last;
						push @cats, $category;
					}
				}
			}
		}
		if ($novel_transcript ==1)
		{
			$TYPE = "-\t-\t-";
			$category = "Novel_transcript";
			push @cats, $category;
		}
		%catz_hash = "";
		if ($protein_coding_seq =~ /$peptide_seq/)
		{
			@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
			foreach $tstars (@twinkle_stars)
			{
				print PRT_OUT "$tstars\t$geno\t$TYPE\t$splice\tCoding\n";
			}
			print PRT_GTF_OUT "$geno\n";
		}
		else
		{
			foreach $cat (@cats)
			{
				$catz_hash{$cat} = 1;
			}
			if (exists $catz_hash{"N_terminal protein extension"}) 
			{ 
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{
					print OUT "$tstars\t$geno\t$TYPE\t$splice\tN_terminal protein extension\n";
				}
				print GTF_OUT "$geno\n";
			}
			elsif (exists $catz_hash{"C_terminal protein extension"})
			{
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{
					print OUT "$tstars\t$geno\t$TYPE\t$splice\tC_terminal protein extension\n";
				}
				print GTF_OUT "$geno\n";
			}
			elsif (exists $catz_hash{"5' UTR"})
			{
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{
					print OUT "$tstars\t$geno\t$TYPE\t$splice\t5' UTR\n";
				}
				print GTF_OUT "$geno\n";
			}
			elsif (exists $catz_hash{"3' UTR"})
			{
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{
					print OUT "$tstars\t$geno\t$TYPE\t$splice\t3' UTR\n";
				}
				print GTF_OUT "$geno\n";
			}
			elsif (exists $catz_hash{"Alternate frame of translation"})
			{
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{	
					print OUT "$tstars\t$geno\t$TYPE\t$splice\tAlternate frame of translation\n";
				}
				print GTF_OUT "$geno\n";
			}
			elsif (exists $catz_hash{"5' Exon extension"})
			{
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{
					print OUT "$tstars\t$geno\t$TYPE\t$splice\t5' Exon extension\n";
				}
				print GTF_OUT "$geno\n";
			}
			elsif (exists $catz_hash{"3' Exon extension"})
			{
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{
					print OUT "$tstars\t$geno\t$TYPE\t$splice\t3' Exon extension\n";
				}
				print GTF_OUT "$geno\n";
			}
			elsif (exists $catz_hash{"Novel Exon"})
			{
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{
					print OUT "$tstars\t$geno\t$TYPE\t$splice\tNovel Exon\n";
				}
				print GTF_OUT "$geno\n";
			}
			elsif (exists $catz_hash{"5' Gene extension"})
			{
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{
					print OUT "$tstars\t$geno\t$TYPE\t$splice\t5' Gene extension\n";
				}
				print GTF_OUT "$geno\n";
			}
			elsif (exists $catz_hash{"3' Gene extension"})
			{
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{
					print OUT "$tstars\t$geno\t$TYPE\t$splice\t3' Gene extension\n";
				}
				print GTF_OUT "$geno\n";
			}
			elsif (exists $catz_hash{"Novel_transcript"})
			{
				@twinkle_stars = split (/\*/, $req_peptide_seq{$peptide_seq});
				foreach $tstars (@twinkle_stars)
				{
					print OUT "$tstars\t$geno\t$TYPE\t$splice\tNovel_transcript\n";
				}
				print GTF_OUT "$geno\n";
			}
		}
	}
}
open (NOT_MAPPED, ">$notmapped_output") or die $!;
while (($ap, $bp) = each %req_peptide_seq)
{
	
	if (!exists $mapped_peptides{$ap})
	{
		@crazy_stars = split (/\*/, $bp);
		$jstars = join ("\t", @crazy_stars);
		@tab_stars = split (/\t/, $crazy_stars[0]);
		print NOT_MAPPED "$jstars\n";
	}
}
