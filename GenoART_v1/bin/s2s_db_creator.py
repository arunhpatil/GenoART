from itertools import groupby
import sys
import os
import re

#function to iter the fasta file
def read_fasta(fasta_name):
    read_file = open(fasta_name)
    faiter = (x[1] for x in groupby(read_file, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            input_file = sys.argv[1]
        ##        idx = input_file.split('.')[-1]
        ##        idx_file = input_file.index(idx)
        ##        file_name_op =  str(input_file[:idx_file - 1]).replace('_pep', '')
        ##            os.rename(input_file, input_file.split('.')[0] + '_' + 'pep.fasta')
            file_name_op = input_file.split('.')[0]
            input_file = input_file.split('.')[0] + '_pep.' + input_file.split('.')[1]
            print "Checking the database for redundancy \n" 
            s2s_file = open(file_name_op + '_s2s.fasta', 'w') # creating the stop to stop sequence fasta file
            m2s_file = open(file_name_op + '_m2s.fasta', 'w') # creating Met to stop sequecne fasta file
            combined_file = open(file_name_op + '_combined.fasta', 'w') # creating the combined(stop to stop and Met to stop sequence) fasta file
            dicts = {}

            #removing the redundancy
            read_fasta_file = read_fasta(input_file)
            for rows in read_fasta_file:
                if rows[1].rstrip() not in dicts:
                    dicts[rows[1].rstrip()] = rows[0].rstrip()

            print "Redundancy check complete!!! \nWriting data to file: \n"

            #iteration and writing output to file
            for k, v in dicts.iteritems():
                if 'M' in k: # checking if 'M' present in dictionary key
                    if len(k[k.index('M'):]) > 6: # checking if length of M to stop has more than 6 amino acid lengths
                        s2s_file.write ('>' + v.rstrip() + '\n' + k + '\n')
                        m2s_file.write ('>' + v.rstrip().replace('|ref', '_m2m|ref') + '\n' + k[k.index('M'):] + '\n')
                        combined_file.write ('>' + v.rstrip() + '\n' + k + '\n')
                        combined_file.write ('>' + v.rstrip().replace('|ref', '_m2m|ref') + '\n' + k[k.index('M'):] + '\n')
                else:
                    s2s_file.write ('>' + v.rstrip() + '\n' + k + '\n')
                    combined_file.write ('>' + v.rstrip() + '\n' + k + '\n')

                #closing the files
            s2s_file.close()        
            m2s_file.close()
            combined_file.close()
        except Exception as e:
            error_log = open('Error.log', 'w')
            print "Please check error log file."
            error_log.write(str(type(e).__name__) + ' Occured at line ' +  str(format(sys.exc_info()[-1].tb_lineno)) + '\n')
    else:
        print "\nUSAGE: s2s_db_creator.exe <input_file.fasta>"
        print "\nExample: s2s_db_creator.exe chr22_pep.fasta\n"
        print "For more details refer to : http://genoart.inhouseprotocols.com/protocols.php\n"
