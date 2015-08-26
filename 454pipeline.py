#!/usr/bin/env python

import os
import re
import sys
import glob
import subprocess
import pipeline_functions2

def main(argv):
    if len(argv) == 5:
        run_pipeline(argv[1], argv[2], argv[3], argv[4])
    elif len(argv) == 6:
        run_pipeline(argv[1], argv[2], argv[3], argv[4], argv[5])
    else:
        print 'usage: python 454pipeline.py folder species trim_threshold overall_quality_threshold primers'
        sys.exit()


def run_pipeline(folder, species, trim_threshold, overall_quality_threshold, primers = False):
    for file_name in glob.glob(folder + '*'):
        pipeline(file_name, species, trim_threshold, overall_quality_threshold, primers)


def pipeline(filename, species, trim_threshold, overall_quality_threshold, primers = False):

    #dissect input filename  
    f = re.compile('(/.*/)[^/]*/([^/]*)\.')
    file_process = f.match(filename)
    folder = file_process.group(1)
    name = file_process.group(2)

    #create directory and filenames for output
    results_dir = folder + name
    cleaned_dir = results_dir + '/cleaned'
    split_dir = cleaned_dir + '/split'
    vidjil_dir = results_dir + '/vidjil'
    clusters_dir = results_dir + '/clusters'
    imgt_dir = results_dir + '/imgt'
    summary_dir = results_dir + '/summary'
    homo_removed = cleaned_dir + '/' + name + '_homo_removed'
    cleaned = cleaned_dir + '/' + name + '_cleaned.fasta'
    collapsed = cleaned_dir + '/' + name + '_collapsed'
    single_removed = cleaned_dir + '/' + name +'_single_removed.fasta'
    imgt_input = imgt_dir + '/' + name

    #make directories for output
    os.system('mkdir ' + results_dir)
    os.system('mkdir ' + cleaned_dir)
    os.system('mkdir ' + vidjil_dir)
    os.system('mkdir ' + imgt_dir)

    if not glob.glob(cleaned):
        #filter reads by length and complexity
        if not glob.glob(homo_removed + '.fastq'):
            os.system('perl /d/as5/s/prinseq/prinseq-lite-0.20.4/prinseq-lite.pl -fastq ' + filename + \
                     ' -lc_method entropy -lc_threshold 70 -min_len 300 -max_len 600 -out_format 3 -line_width 0 -out_bad null -out_good ' + homo_removed + ' &> ' + cleaned_dir + '/prinseq.txt')
            print 'Reads filtered by length and low complexity reads removed'

        #remove primers, trim by quality and produce QC graphs
        if primers:
            if not glob.glob(cleaned_dir + '/QC.unpaired.trimmed.fastq'):
                os.system('perl /d/as5/u/la001/software/FaQCs.pl -u ' + homo_removed + '.fastq -d ' + cleaned_dir + ' -artifactFile ' + folder + '/primers/primers.fasta -5trim_off -q '\
                          + trim_threshold)
                print 'Primers removed, 3\' ends trimmed by quality, QC graphs produced'
        elif not glob.glob(cleaned_dir + '/QC.unpaired.trimmed.fastq'):
            os.system('perl /d/as5/u/la001/software/FaQCs.pl -u ' + homo_removed + '.fastq -d ' + cleaned_dir + ' -5trim_off -q 20 -min_L 300')
            print '3\' ends trimmed by quality, QC graphs produced'

       #filter reads by quality
        if not glob.glob(cleaned):
            subprocess.call('cat ' + cleaned_dir + '/QC.unpaired.trimmed.fastq | split -l 40000000 -a 2 - ' + split_dir, shell = True)
            for f in glob.glob(split_dir + 'a*'):
               subprocess.call('/d/user5/wlees01/bin/usearch -threads 15 -fastq_filter ' + f + ' -fastq_maxee ' + overall_quality_threshold + ' -fastaout ' + f + '.fasta &> ' + cleaned_dir + '/usearch.txt' , shell = True)
            subprocess.call('cat ' + split_dir + 'a*.fasta > ' + cleaned, shell = True)
            print 'Reads filtered by quality'

        #remove temporary files
        os.system('rm ' + homo_removed + '.fastq')
        os.system('rm '+ cleaned_dir + '/QC.unpaired.trimmed.fastq')
        print 'Tempory files removed'
             
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

 
if __name__=="__main__":
    main(sys.argv)                                    
        



        






