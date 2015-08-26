#Functions to determine statistics associated with the isoelectric point of CDR3s with a certain property. 
#A textfile listing IMGT output folders and output file name are taken as input for 'run_pi'. 
#Additionally a fuction to determine whether CDR3s satisfy a certain property must be supplied as an input argument.
#These are included in helper_functions.py.
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID

import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import helper_functions
import fractions


def find_pi_stats(folders, outfile, condition, regex = None, frequency = False):
    mean_list = []
    for folder in folders:
        with open(folder[0] + '/6_Junction.txt') as f:
            pi_all = 0
            total_seqs = 0
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'] == 'productive' and condition(row['CDR3-IMGT (AA)']):
                    try:
                        pi = float(row['pI'])
                        if frequency:
                            pat = re.compile(regex)
                            info = pat.match(row['Sequence ID'])
                            freq = int(info.group(1))
                        else:
                            freq = 1
                        total_seqs += freq
                        pi_all += pi * freq
                    except:
                        pass
            mean_list.append(pi_all/float(total_seqs))
    with open(outfile + '_means.txt', 'w') as out:
        for item in mean_list:
            out.write(str(item) +'\n')
    with open(outfile + '.txt', 'w') as out:
        out.write('mean CDR3 pi,standard deviation\n')
        out.write(str(np.mean(mean_list)) + ',' + str(np.std(mean_list)))

def run_pi(infile, outfile, condition, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    find_pi_stats(folders, outfile, condition, regex, frequency)


