from itertools import groupby
from itertools import islice
from string import maketrans
import re
import sys
import copy

print "##################################################"
print "### GenoART: Genome Annotation Refinement Tool ###"
print "###     GenoART v1.0.0 2017, \"Arun H Patil\"    ###"
print "##################################################"
# GenoART is free software: you can redistribute it and/or modify
# it and, is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

def read_fasta(fasta_name):
    read_file = open(fasta_name)
    faiter = (x[1] for x in groupby(read_file, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

              
complementrary_input = "ATGC"
complemtary_output = "TACG"
translated_seqs = maketrans(complementrary_input, complemtary_output)
dict_conv = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
#dicts_gtf = {}
dicts_exon_gtf = {}
dicts_CDS_gtf = {}
dicts_mrna_gtf = {}
def load_gtf(infile_gff):
    flag_exon = 0
    flag_cds = 0
    flag_mrna = 0 
    for each_line in open(infile_gff):
            split_each_line = each_line.split('\t')
            if not each_line.startswith('#'):
                if split_each_line[2] == 'exon':
                    flag_exon = 1
                    if re.search('transcript_id (.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('transcript_id (.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_exon_gtf:
                            dicts_exon_gtf[chromosome] = {}
                            #dicts_exon_gtf[chromosome] = {}
                        if transcript_id not in dicts_exon_gtf[chromosome]:
                            dicts_exon_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_exon_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                            
                    elif re.search('ID=(.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('ID=(.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_exon_gtf:
                            dicts_exon_gtf[chromosome] = {}
                        if transcript_id not in dicts_exon_gtf[chromosome]:
                            dicts_exon_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_exon_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                            
                    elif re.search('ID=exon:(.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('ID=exon:(.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_exon_gtf:
                            dicts_exon_gtf[chromosome] = {}
                        if transcript_id not in dicts_exon_gtf[chromosome]:
                            dicts_exon_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_exon_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                            
                    elif re.search('ID=exon_(.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('ID=exon_(.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_exon_gtf:
                            dicts_exon_gtf[chromosome] = {}
                        if transcript_id not in dicts_exon_gtf[chromosome]:
                            dicts_exon_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_exon_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)

                    else:
                        tra_id = 1
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], 'transcript_id_' + str(tra_id), split_each_line[6]
                        if chromosome not in dicts_exon_gtf:
                            dicts_exon_gtf[chromosome] = {}
                        if transcript_id not in dicts_exon_gtf[chromosome]:
                            dicts_exon_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_exon_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                        tra_id+=1
                    

                if split_each_line[2] == 'CDS':
                    flag_cds = 2
                    if re.search('transcript_id (.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('transcript_id (.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_CDS_gtf:
                            dicts_CDS_gtf[chromosome] = {}
                            #dicts_exon_gtf[chromosome] = {}
                        if transcript_id not in dicts_CDS_gtf[chromosome]:
                            dicts_CDS_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_CDS_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                            
                    elif re.search('ID=(.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('ID=(.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_CDS_gtf:
                            dicts_CDS_gtf[chromosome] = {}
                        if transcript_id not in dicts_CDS_gtf[chromosome]:
                            dicts_CDS_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_CDS_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                            
                    elif re.search('ID=exon:(.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('ID=exon:(.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_CDS_gtf:
                            dicts_CDS_gtf[chromosome] = {}
                        if transcript_id not in dicts_CDS_gtf[chromosome]:
                            dicts_CDS_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_CDS_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                            
                    elif re.search('ID=exon_(.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('ID=exon_(.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_CDS_gtf:
                            dicts_CDS_gtf[chromosome] = {}
                        if transcript_id not in dicts_CDS_gtf[chromosome]:
                            dicts_CDS_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_CDS_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
    
                    else:
                        tra_id = 1
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], 'transcript_id_' + str(tra_id), split_each_line[6]
                        if chromosome not in dicts_CDS_gtf:
                            dicts_CDS_gtf[chromosome] = {}
                        if transcript_id not in dicts_CDS_gtf[chromosome]:
                            dicts_CDS_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_CDS_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                        tra_id+=1
                    
                    
                if split_each_line[2] == 'mRNA':
                    flag_mrna = 3
                    if re.search('transcript_id (.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('transcript_id (.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_mrna_gtf:
                            dicts_mrna_gtf[chromosome] = {}
                            #dicts_exon_gtf[chromosome] = {}
                        if transcript_id not in dicts_mrna_gtf[chromosome]:
                            dicts_mrna_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_mrna_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                            
                    elif re.search('ID=(.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('ID=(.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_mrna_gtf:
                            dicts_mrna_gtf[chromosome] = {}
                        if transcript_id not in dicts_mrna_gtf[chromosome]:
                            dicts_mrna_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_mrna_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                            
                    elif re.search('ID=exon:(.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('ID=exon:(.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_mrna_gtf:
                            dicts_mrna_gtf[chromosome] = {}
                        if transcript_id not in dicts_mrna_gtf[chromosome]:
                            dicts_mrna_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_mrna_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                            
                    elif re.search('ID=exon_(.*?)\;', split_each_line[8]):
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], split_each_line[8], split_each_line[6]
                        transcript_id = re.search('ID=exon_(.*?)\;', desc_id).group(1)
                        if chromosome not in dicts_mrna_gtf:
                            dicts_mrna_gtf[chromosome] = {}
                        if transcript_id not in dicts_mrna_gtf[chromosome]:
                            dicts_mrna_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_mrna_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)

                    else:
                        tra_id = 1
                        chromosome, coord_start, coord_end, desc_id, strand = split_each_line[0], split_each_line[3], split_each_line[4], 'transcript_id_' + str(tra_id), split_each_line[6]
                        if chromosome not in dicts_mrna_gtf:
                            dicts_mrna_gtf[chromosome] = {}
                        if transcript_id not in dicts_mrna_gtf[chromosome]:
                            dicts_mrna_gtf[chromosome][transcript_id] = [coord_start + ";" + coord_end + ";" +strand]
                        else:
                            dicts_mrna_gtf[chromosome][transcript_id].append(coord_start + ";" + coord_end + ";" +strand)
                        tra_id+=1
    return flag_exon, flag_cds, flag_mrna
                    
def no_gtf(infile):
    cnt = 100001
    idx = infile.split('.')[-1]
    idx_file = infile.index(idx)
    file_name_op =  infile[:idx_file - 1]
    print "\nInput: \n\tGenome file is: "+ sys.argv[1] +"\n\tGTF file is: NULL\n\nFetching sequence is in progress!!!\n"
    write_file = open(file_name_op + '_op.fasta', 'w')
    read_file = read_fasta(infile)
    for rows in read_file:
        write_file.write('>' + str(cnt) + ' ' + 'Trans_' + str(cnt) + ' ' + rows[0].split(' ')[0] + '+' +' ' + '1-' + str(len(rows[1])) + '\n')
        write_file.write(rows[1].upper() + '\n')
        cnt+=1
        write_file.write('>' + str(cnt) + ' ' + 'Trans_' + str(cnt)+ ' ' + rows[0].split(' ')[0] + '-' +' ' + '1-' + str(len(rows[1])) + '\n')
        reverse_sequence_d = rows[1][::-1]
        reverse_sequence = reverse_sequence_d.upper()
        complementary_seqs = reverse_sequence.translate(translated_seqs)
        write_file.write(complementary_seqs + '\n')
        cnt+=1

dicts_gtf = {}
if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            chr_s = ''
            for i in open(sys.argv[1]):
                if i.startswith('>'):
                    ids = i.split(' ')
                    if len(ids) != 0:
                        chr_s = chr_s + ids[0].rstrip()[1:] + '; '
                    else:
                        chr_s = chr_s + i.rstrip()[1:] + '; '
            input_gtf = sys.argv[2]

            if input_gtf != '0':
                a = load_gtf(input_gtf)  # Input GTF file
                print "\nInput: \n\tGenome file is: "+ sys.argv[1] +"\n\tGTF file is: "+input_gtf+"\n\nFetching sequence is in progress!!!\n\nThe annotation file contains the following chromosomes/scaffolds:"
                print chr_s + '\n'
                load_fasta_file = sys.argv[1]
                idx = load_fasta_file.split('.')[-1]
                idx_file = load_fasta_file.index(idx)
                file_name_op =  load_fasta_file[:idx_file - 1]
                write_file1 = open(file_name_op + '_op.fasta', 'w')
                cnt = 1000 #chr id starts from 1000
                store_seq_negative = ''
                store_seq_positive = ''
                store_seq_positive_d = ''
                store_seq_negative_d = ''
                store_seq_unknown = ''
                store_seq_unknown_d = ''
                if a[0] == 1:
                    dicts_gtf = copy.deepcopy(dicts_exon_gtf)
                elif a[1] == 2:
                    dicts_gtf = copy.deepcopy(dicts_CDS_gtf)
                elif a[2] == 3:
                    dicts_gtf = copy.deepcopy(dicts_mrna_gtf)
                for rows in read_fasta(load_fasta_file): # Input Fasta genome sequence
                    if '|' in rows[0]:
                        rows_pipe = rows[0].split('|')[0].rstrip()
                    else:
                        rows_pipe = rows[0].split(' ')[0].rstrip()
                    if rows_pipe in dicts_gtf: # rows[0]
                        for key, val in dicts_gtf[rows_pipe].iteritems(): # rows[0]
                            if '-' in str(val):
                                st_cd = int(val[0].split(';')[0]) # frist item of the list sepeated by ";'
                                ed_cd = int(val[-1].split(';')[0]) # last item of the list seperated by ";"
                                if st_cd < ed_cd:    
                                    for iter_each_item in val:
                                        split_coords = iter_each_item.split(';')
                                        coord_start_pos, coord_end_pos = int(split_coords[0]), int(split_coords[1])
                                        store_seq_negative_d = store_seq_negative_d + rows[1][coord_start_pos -1: coord_end_pos]
                                        store_seq_negative = store_seq_negative_d.upper()
                                    reversed_coordinates = [x for x in val[::-1]]
                                    reverse_sequence = store_seq_negative[::-1]
                                    complementary_seqs = reverse_sequence.translate(translated_seqs)
                                    write_file1.write('>' + str(cnt) + ' ' + key.replace('"', '') + ' ' + rows_pipe + '-' + ' ' + str(reversed_coordinates).replace(';-', '').replace("['", '').replace("']", '').replace("'", '').replace(";", "-").replace(' ', '') + '\n' + complementary_seqs + '\n')
                                    store_seq_negative = ''
                                    store_seq_negative_d = ''
                                    cnt+=1
                                else:
                                    for iter_each_item in val[::-1]:
                                        split_coords = iter_each_item.split(';')
                                        coord_start_pos, coord_end_pos = int(split_coords[0]), int(split_coords[1])
                                        store_seq_negative_d = store_seq_negative_d + rows[1][coord_start_pos -1: coord_end_pos]
                                        store_seq_negative = store_seq_negative_d.upper()
                                    reversed_coordinates = [x for x in val[::-1]]
                                    reverse_sequence = store_seq_negative[::-1]
                                    complementary_seqs = reverse_sequence.translate(translated_seqs)
                                    write_file1.write('>' + str(cnt) + ' ' + key.replace('"', '') + ' ' + rows_pipe + '-' + ' ' + str(reversed_coordinates).replace(';-', '').replace("['", '').replace("']", '').replace("'", '').replace(";", "-").replace(' ', '') + '\n' + complementary_seqs + '\n')
                                    store_seq_negative = ''
                                    store_seq_negative_d = ''
                                    cnt+=1
                                    
                            elif '+' in str(val):
                                for iter_each_item1 in val:
                                    split_coords1 = iter_each_item1.split(';')
                                    coord_start_pos1, coord_end_pos1 = int(split_coords1[0]), int(split_coords1[1])
                                    store_seq_positive_d = store_seq_positive_d + rows[1][coord_start_pos1 - 1: coord_end_pos1]
                                    store_seq_positive = store_seq_positive_d.upper()
                                write_file1.write('>' + str(cnt) + ' '  + key.replace('"', '') + ' '  + rows_pipe + '+' + ' ' + str(val).replace(';+', '').replace("['", '').replace("']", '').replace("'", '').replace(";", "-").replace(' ', '')+ '\n' + store_seq_positive + '\n')
                                store_seq_positive = ''
                                store_seq_positive_d = ''
                                cnt+=1
                                      
                            else:
                                st_cd2 = int(val[0].split(';')[0]) # frist item of the list sepeated by ";'
                                ed_cd2 = int(val[-1].split(';')[0]) # last item of the list seperated by ";"
                                if st_cd2 < ed_cd2:
                                    for iter_each_item2 in val:
                                        split_coords2 = iter_each_item2.split(';')
                                        coord_start_pos2, coord_end_pos2 = int(split_coords2[0]), int(split_coords2[1])
                                        store_seq_unknown_d = store_seq_unknown_d + rows[1][coord_start_pos2 - 1: coord_end_pos2]
                                        store_seq_unknown = store_seq_unknown_d.upper()
                                    reversed_coordinates1 = [x for x in val[::-1]]
                                    reverse_sequence1 = store_seq_unknown[::-1]
                                    complementary_seqs1 = reverse_sequence1.translate(translated_seqs)
                                    write_file1.write('>' + str(cnt) + ' '  + key.replace('"', '') + ' '  + rows_pipe + '+' + ' ' + str(val).replace(';.', '').replace("['", '').replace("']", '').replace("'", '').replace(";", "-").replace(' ', '')+ '\n' + store_seq_unknown + '\n')
                                    cnt+=1
                                    write_file1.write('>' + str(cnt) + ' '  + key.replace('"', '') + ' '  + rows_pipe + '-' + ' ' + str(val).replace(';.', '').replace("['", '').replace("']", '').replace("'", '').replace(";", "-").replace(' ', '')+ '\n' + complementary_seqs1 + '\n')
                                    store_seq_unknown = ''
                                    store_seq_unknown_d = ''
                                    cnt+=1
                                else:
                                    for iter_each_item2 in val[::-1]:
                                        split_coords2 = iter_each_item2.split(';')
                                        coord_start_pos2, coord_end_pos2 = int(split_coords2[0]), int(split_coords2[1])
                                        store_seq_unknown_d = store_seq_unknown_d + rows[1][coord_start_pos2 - 1: coord_end_pos2]
                                        store_seq_unknown = store_seq_unknown_d.upper()
                                    reversed_coordinates1 = [x for x in val[::-1]]
                                    reverse_sequence1 = store_seq_unknown[::-1]
                                    complementary_seqs1 = reverse_sequence1.translate(translated_seqs)
                                    write_file1.write('>' + str(cnt) + ' '  + key.replace('"', '') + ' '  + rows_pipe + '+' + ' ' + str(val).replace(';.', '').replace("['", '').replace("']", '').replace("'", '').replace(";", "-").replace(' ', '')+ '\n' + store_seq_unknown + '\n')
                                    cnt+=1
                                    write_file1.write('>' + str(cnt) + ' '  + key.replace('"', '') + ' '  + rows_pipe + '-' + ' ' + str(val).replace(';.', '').replace("['", '').replace("']", '').replace("'", '').replace(";", "-").replace(' ', '')+ '\n' + complementary_seqs1 + '\n')
                                    store_seq_unknown = ''
                                    store_seq_unknown_d = ''
                                    cnt+=1
                                                   
                    print "Sequence fetched for: " + rows_pipe # Prints the number of chromosomes processed
                write_file1.close()
                print "\nFetching sequence is complete!!!\n"
            else:
                run_scr = no_gtf(sys.argv[1])
                print "\nFetching sequence is complete!!!\n"
        except Exception as e:
            error_log = open('Error.log', 'w')
            ##        error_log.write('Error on line ' + str(format(sys.exc_info()[-1].tb_lineno)) + '\n')
            ##        error_log.write(str(e) + '\n')
            ##        if format(sys.exc_info()[-1].tb_lineno) == '65':
            ##            error_log.write(str(type(e).__name__) + ' Occured. ' + "Input GTF and Fasta files are missing.")
            ##        else:
            error_log.write(str(type(e).__name__) + ' Occured at line ' +  str(format(sys.exc_info()[-1].tb_lineno)) + '\n')
    else:
         print "\nUSAGE: fetch_seq.exe <Genome_file.fasta> <Reference_GTF.gtf>"
         print "\nExample: fetch_seq.exe chr22.fasta chr_hs.gtf\n"
         print "For more details refer to : http://genoart.inhouseprotocols.com/protocols.php\n"
         print "Please check Error.log file \n"
                


        

    
                
