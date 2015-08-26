#Script to analyse CDR3 lengths. 
#A textfile listing IMGT output folders and output file name are taken as input. 
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID.

import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import helper_functions

def main(argv):
    if len(argv) == 3:
        run_lengths(argv[1], argv[2])
    elif len(argv) == 5:
	run_lengths(argv[1], argv[2], argv[3], argv[4])
    else:
        print 'usage: CDR3_length.py IMGT_folder_list outfile '
	print 'optional arguments: frequency_regex frequency(TRUE/FALSE)'
        sys.exit()



def find_length_stats(folders, outfile, regex = None, frequency = False):
    mean_list = []
    for folder in folders:
        with open(folder[0] + '/5_AA-sequences.txt') as f:
            lengths = 0
            total_seqs = 0
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'] == 'productive':
                    length = len(row['CDR3-IMGT'])
                    if frequency:
                        pat = re.compile(regex)
                        info = pat.match(row['Sequence ID'])
                        freq = int(info.group(1))
                    else:
                        freq = 1
                    total_seqs += freq
                    lengths += length * freq
            mean_list.append(lengths/float(total_seqs))
    with open(outfile + '_means.txt', 'w') as out:
        for item in mean_list:
            out.write(str(item) + ',')
    with open(outfile, 'w') as out:
        out.write('mean CDR3 length,standard deviation\n')
        out.write(str(np.mean(mean_list)) + ',' + str(np.std(mean_list)))

def run_lengths(infile, outfile, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    find_length_stats(folders, outfile, regex, frequency)

if __name__=="__main__":
     main(sys.argv) 



