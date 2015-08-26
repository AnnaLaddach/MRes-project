import re
import glob
import sys
import os
import move_barcode

def main(argv):
    if len(argv) == 2:
        run_illumina_pipeline(argv[1])
    else:
        print 'usage: flu_pipeline.py folder'
        sys.exit()


def run_illumina_pipeline(folder):
    for file_name in glob.glob(folder + '*_1.fastq'):
        illumina_pipeline(file_name)


def illumina_pipeline(infile):

    pat = re.compile('(^.+)/.*?/(.*?)\.fastq')
    info = pat.search(infile)
    folder = info.group(1)
    name = info.group(2)

    results_dir = folder + '/' + name
    cleaned_dir = results_dir + '/cleaned'
    
    os.system('mkdir ' + results_dir)
    os.system('mkdir ' + cleaned_dir)

    if not glob.glob(results_dir + '/' + name[:-2] + '.fasta.cdr'):

	#move barcode to forward read
        move_barcode.move_barcode(infile, infile[:-7] + '2.fastq', cleaned_dir + '/')
        
        #filter reads by quality
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/FilterSeq.py trimqual -s ' + cleaned_dir + '/R1.fastq -q 20 --nproc 4 --outname R1 --clean --outdir ' + cleaned_dir + ' &> ' + cleaned_dir + '/quality_R1.txt' ) 
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/FilterSeq.py trimqual -s ' + cleaned_dir + '/R2.fastq -q 20 --nproc 4 --outname R2 --outdir ' + cleaned_dir + ' &> ' + cleaned_dir + '/quality_R2.txt')
        print 'reads cleaned'

        #mask primers and identify UIDs 
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/MaskPrimers.py score -s ' + cleaned_dir + '/R1_trimqual-pass.fastq -p /l_mnt/ssd1/u/la001/primers/forward.fasta --mode cut --start 16 --barcode --maxerror 0.2 --nproc 4 --clean &> ' + cleaned_dir + '/R1_mask_primers.txt')
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/MaskPrimers.py score -s ' + cleaned_dir + '/R2_trimqual-pass.fastq -p /l_mnt/ssd1/u/la001/primers/reverse.fasta --mode cut --start 0 --maxerror 0.2 --nproc 4 --clean &> ' + cleaned_dir + '/R2_mask_primers.txt')
        print 'primers masked, UIDs identified'

        #assign UIDs to R2 reads
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/PairSeq.py -1 ' + cleaned_dir + '/R1*primers-pass.fastq -2 ' + cleaned_dir + '/R2*primers-pass.fastq -f BARCODE --coord sra --clean &> ' + cleaned_dir + '/pair_seqs.txt')
        print 'UIDs assigned'

        #align sequence start positions by primer alignments
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/AlignSets.py offset -s' + cleaned_dir + '/R1_trimqual-pass_primers-pass_pair-pass.fastq -d /l_mnt/ssd1/u/la001/primers/forward_offsets.tab --bf BARCODE --pf PRIMER --nproc 4 &> ' + cleaned_dir + '/R1_seqs_aligned.txt')
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/AlignSets.py offset -s' + cleaned_dir + '/R2_trimqual-pass_primers-pass_pair-pass.fastq -d /l_mnt/ssd1/u/la001/primers/reverse_offsets.tab --bf BARCODE --pf PRIMER --nproc 4 &> ' + cleaned_dir + '/R2_seqs_aligned.txt')
        print 'sequences aligned'

        #build consensus seqs
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/BuildConsensus.py -s ' + cleaned_dir + '/R1*align-pass.fastq --bf BARCODE --pf PRIMER --prcons 0.6 -q 20 --maxdiv 0.1 -n 2 --nproc 4 &> ' + cleaned_dir + '/consensus_R1.txt')
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/BuildConsensus.py -s ' + cleaned_dir + '/R2*align-pass.fastq --bf BARCODE --pf PRIMER --prcons 0.6 -q 20 --maxdiv 0.1 -n 2 --nproc 4 &> ' + cleaned_dir + '/consensus_R2.txt')
        print 'consensi built'

        #assemble paired-end reads
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/AssemblePairs.py align -1 ' + cleaned_dir + '/R1*consensus-pass.fastq -2 ' + cleaned_dir + '/R2*consensus-pass.fastq --1f CONSCOUNT --2f PRCONS CONSCOUNT --coord presto --rc tail --maxerror 0.05 --minlen 10 --nproc 4 '\
                  + '--outname Assembled --clean &> ' + cleaned_dir + '/assembled.txt')
        print 'sequences assembled'

        #remove sequences with gaps
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/FilterSeq.py missing -s ' + cleaned_dir + '/Assembled_assemble-pass.fastq -n 0 --inner &> ' + cleaned_dir + '/missing_filtered.txt')
        print 'gapped sequences removed'
      
        #collapse reads
        os.system('/d/as2/u/wlees01/fastx/bin/fastx_collapser -i ' + cleaned_dir + '/Assembled*missing-pass.fastq -o ' + results_dir + '/' + name[:-2] + '.fasta -Q 33')
        print 'reads collapsed'

        #find CDR3s
        os.system('python /d/as5/u/la001/software/AbMiningToolbox/cdr3_search2.py ' + results_dir + '/' + name[:-2] + '.fasta')
        print 'CDR3s found'

	#remove uneeded files
    if glob.glob(results_dir + '/' + name[:-2] + '.fasta'):
        os.system('rm -r ' + cleaned_dir + '/*.fastq')
        print 'fastq files removed'

        

if __name__=="__main__":
     main(sys.argv)   

     

        

        

        

        
        
        
        
    
