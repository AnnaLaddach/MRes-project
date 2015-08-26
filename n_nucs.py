#Function to determine statistics associated with the n-nucleotide usage of CDR3s with a certain property. 
#A textfile listing IMGT output folders and output file name are taken as input for 'run_n_nucs'. 
#Additionally a fuction to determine whether CDR3s satisfy a certain property must be supplied as an input argument.
#These are included in helper_functions.py.
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID

import csv
import re
import math
import numpy as np
import matplotlib.pyplot as plt
import helper_functions
import fractions


def find_n_nucs_stats(folders, outfile, condition, regex = None, frequency = False):
    mean_list = []
    for folder in folders:
        with open(folder[0] + '/6_Junction.txt') as f:
            n_all = 0
            total_seqs = 0
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                try:
                    if row['Functionality'] == 'productive' and condition(row['CDR3-IMGT (AA)']):
                        n_nucs = int(row['N1-REGION-nt nb']) + int(row['N2-REGION-nt nb'])\
                                 + int(row['N3-REGION-nt nb']) + int(row['N4-REGION-nt nb'])\
                                 + int(row['N-REGION-nt nb'])
                        if frequency:
                            pat = re.compile(regex)
                            info = pat.match(row['Sequence ID'])
                            freq = int(info.group(1))
                        else:
                            freq = 1
                        total_seqs += freq
                        n_all += n_nucs * freq
                except:
                    pass
            try:    
                mean_list.append(n_all/float(total_seqs))
                print mean_list
            except:
                pass
    with open(outfile + '_means.txt', 'w') as out:
        for item in mean_list:
            out.write(str(item) +'\n')
    with open(outfile + '.txt', 'w') as out:
        out.write('mean CDR3 n nucleotides,standard deviation\n')
        out.write(str(np.mean(mean_list)) + ',' + str(np.std(mean_list)))

def run_n_nucs(infile, outfile, condition, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    find_n_nucs_stats(folders, outfile, condition, regex, frequency)


