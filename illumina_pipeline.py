import os
import re
import split_reads2
import sys
import glob
import pipeline_functions2
import subprocess

def main(argv):
    if len(argv) == 4:
        run_illumina_pipeline(argv[1], argv[2], argv[3])
    else:
        print 'usage: python illumina_pipeline_3.py folder species overall_quality_threshold'
        sys.exit()


def run_illumina_pipeline(folder, species, overall_quality_threshold):
    for file_name in glob.glob(folder + '*_1.fastq'):
        illumina_pipeline(file_name, species, overall_quality_threshold)
    

def illumina_pipeline(infile, species, overall_quality_threshold):
    pat = re.compile('(^.+)/.*?/(.*?)\.fastq')
    info = pat.search(infile)
    folder = info.group(1)
    name = info.group(2)
    
    results_dir = folder + '/' + name
    cleaned_dir = results_dir + '/cleaned'
    vidjil_dir = results_dir + '/vidjil'
    imgt_dir = results_dir + '/imgt'
    cleaned = cleaned_dir + '/' + name + '_cleaned.fasta'
    imgt_input = imgt_dir + '/' + name
    split_dir = cleaned_dir + '/split'

    #make output directory
    os.system('mkdir ' + results_dir)
    os.system('mkdir ' + cleaned_dir)
    os.system('mkdir ' + vidjil_dir)
    os.system('mkdir ' + imgt_dir)
    if not glob.glob(cleaned):
        #trim reads by quality
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/FilterSeq.py trimqual -s ' + infile + ' -q 20 --nproc 4 --outname R1 --clean --outdir ' + cleaned_dir + ' &> ' + cleaned_dir + '/quality_R1.txt') 
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/FilterSeq.py trimqual -s ' + infile[:-7] + '2.fastq -q 20 --nproc 4 --outname R2 --outdir ' + cleaned_dir + ' &> ' + cleaned_dir + '/quality_R2.txt') 
	print 'reads trimmed'
        
        #mask primers
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/MaskPrimers.py align -s ' + cleaned_dir + '/R1_trimqual-pass.fastq -p /d/as7/scratch/u/la001/SRP057017/artifacts/forward.fasta --mode cut --maxerror 0.2 --nproc 4 --log PrimerLogR1.log --clean &> ' + cleaned_dir + '/R1_mask_primers.txt')
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/MaskPrimers.py align -s ' + cleaned_dir + '/R2_trimqual-pass.fastq -p /d/as7/scratch/u/la001/SRP057017/artifacts/reverse.fasta --mode cut --maxerror 0.2 --nproc 4 --log PrimerLogR2.log --clean &> ' + cleaned_dir + '/R2_mask_primers.txt') 
	print 'primers removed'	

        #join reads
        os.system('python2.7 /d/as5/u/la001/software/pRESTO_v0.4/AssemblePairs.py align -1 ' + cleaned_dir + '/R1*primers-pass.fastq -2 ' + cleaned_dir + '/R2*primers-pass.fastq --coord sra --rc tail --maxerror 0.05 --minlen 10 --nproc 4 --log AssembleLog.log --outname Assembled --clean &> ' + cleaned_dir + '/assembled.txt')
	print 'reads joined'

        if glob.glob(cleaned_dir + '/Assembled_assemble-pass.fastq') and not glob.glob(cleaned_dir + '/QC.unpaired.trimmed.fastq'):
                os.system('perl /d/as5/u/la001/software/FaQCs.pl -u ' + cleaned_dir + '/Assembled_assemble-pass.fastq' + ' -d ' + cleaned_dir + ' -qc_only') 
		print 'quality control graphs made'      
   
        if glob.glob(cleaned_dir + '/Assembled_assemble-pass.fastq') and not glob.glob(cleaned):
            subprocess.call('cat ' + cleaned_dir + '/Assembled_assemble-pass.fastq | split -l 40000000 -a 2 - ' + split_dir, shell = True)
            for f in glob.glob(split_dir + 'a*'):
               subprocess.call('/d/user5/wlees01/bin/usearch -threads 15 -fastq_filter ' + f + ' -fastq_maxee ' + overall_quality_threshold + ' -fastaout ' + f + '.fasta &> ' + cleaned_dir + '/usearch.txt' , shell = True)
            subprocess.call('cat ' + split_dir + 'a*.fasta > ' + cleaned, shell = True)
            print 'Reads filtered by quality'

    #gather into clones
    if not glob.glob(vidjil_dir + '/output.log'):
        os.system('/d/as2/u/wlees01/vidjil/vidjil -c clones -o ' + vidjil_dir + \
                  ' -G /d/as2/u/wlees01/vidjil/germline/' + species + '/IGH -y all -z 0 -r 2 -d -w 60 ' + \
                  cleaned + ' > ' + vidjil_dir + '/output.log')
        print 'Reads gathered into clones'

    #split fasta output from vidjil into files processable by IGMT, determine number of files.
    if not glob.glob(imgt_input+'*'):
        num_files = pipeline_functions2.split_vidjil_output(vidjil_dir + '/' + name + '_cleaned.vdj.fa', imgt_input)
        print 'IMGT input files prepared'

    #extract CDR3s
    if glob.glob(imgt_input + '*') and not glob.glob(imgt_input + '1.fasta.cdr'):
        for seq_file in glob.glob(imgt_input + '*'):
            os.system('python /d/as5/u/la001/software/AbMiningToolbox/cdr3_search.py ' + seq_file)

    #remove fastq files
    if glob.glob(imgt_input + '*'):
        os.system('rm -r ' + cleaned_dir + '/*.fastq')
      
   
if __name__=="__main__":
    main(sys.argv)  
