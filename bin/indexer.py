import sys
from pyfaidx import Faidx

if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            fa = Faidx(sys.argv[1])
            fa.index
        except Exception as e:
            error_log = open('Error.log', 'w')
            error_log.write(str(type(e).__name__) + ' Occured at line ' +  str(format(sys.exc_info()[-1].tb_lineno)) + '\n')
    else:
         print "\nUSAGE: indexer.exe <Genome_file.fasta>"
         print "\nExample: indexer.exe chr22.fasta\n"
         print "For more details refer to : http://genoart.inhouseprotocols.com/protocols.php\n"
         print "Please check Error.log file \n"
                
        
