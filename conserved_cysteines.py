#Script to analyse and plot the relative frequency of cysteine conservation from IGHD2 germline gene segments in their second reading frame in relation to CDR3 length. Only those D segments assigned enough nucleotides to ensure both cysteines were included in the initial recombination event are used for analysis (>24 nucleotides for IGHD2-2, 2-15 and 2-8; >21 nucleotides for IGHD2-21). A textfile listing IMGT output folders and output file name are taken as input. Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID.

import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import helper_functions
import fractions

class Cdr:
    def __init__(self, length, cysteines, frequency):
        self._length = length
        self._number = 0
        self._cysteines = 0
        self._fraction_conserved = 0
        self._increment(cysteines, frequency)
        
    def _increment(self, cysteines, frequency):
        self._number += frequency
        self._cysteines += cysteines

    def _calculate_stats(self):
        self._fraction_conserved = self._cysteines / float(self._number)

class Cdr_all_files:
    def __init__(self, length, cysteines, number):
        self._length = length
        self._cysteines = []
        self._number =0
        self._files = 0
        self._increment(cysteines, number)
      
    def _increment(self, cysteines, number):
        self._cysteines.append(cysteines)
        self._files += 1
        self._number += number

    def _calculate_stats(self):
        self._cysteines = np.mean(self._cysteines) 

def main(argv):
    if len(argv) == 3:
        plot_conserved_cysteines(argv[1], argv[2])
    elif len(argv) == 5:
	plot_conserved_cysteines(argv[1], argv[2], argv[3], argv[4])
    else:
        print 'usage: conserved_cysteines.py IMGT_folder_list outfile '
	print 'optional arguments: frequency_regex frequency(TRUE/FALSE)'
        sys.exit()


def conserved_cysteines(folders, outfile, regex = None, frequency = False):
    cdr_list_all = []
    for folder in folders:
        with open(folder[0] + '/6_Junction.txt') as f:
            cdr_list = []
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'].find('productive') == 0 and row['D-REGION reading frame'] == '2':
                    try:
                        pat = re.compile('IGHD[^\*]+')
                        info = pat.search(row['D-GENE and allele'])
                        dgene = info.group(0)
                    except:
                        dgene = ''
                    if frequency:
                        pat = re.compile(regex)
                        info = pat.match(row['Sequence ID'])
                        freq = int(info.group(1))
                    else:
                        freq = 1
                    cysteines = 0
                    d2 = False
                    if dgene in ['IGHD2-2','IGHD2-8','IGHD2-15'] and int(row['D-REGION-nt nb']) > 24:
                        d2 = True
                        if helper_functions.spacing_4(row['CDR3-IMGT (AA)']):
                            cysteines = freq
                    if dgene == 'IGHD2-21' and int(row['D-REGION-nt nb']) > 21:
                        d2 = True
                        if helper_functions.spacing_3(row['CDR3-IMGT (AA)']):
                            cysteines = freq
                    if d2:
                        length = len(row['CDR3-IMGT (AA)'])
                        notfound = True
                        for item in cdr_list:
                            if length == item._length:
                                item._increment(cysteines, freq)
                                notfound = False
                        if notfound:
                            cdr_list.append(Cdr(length, cysteines, freq))
        for item in cdr_list:
            item._calculate_stats()
            notfound2 = True
            for item2 in cdr_list_all:
                if item._length == item2._length:
                    item2._increment(item._fraction_conserved, item._number)
                    notfound2 = False
            if notfound2:
                cdr_list_all.append(Cdr_all_files(item._length, item._fraction_conserved, item._number))
    cdr_list_all.sort(key = lambda x: x._length)
    data = []
    with open(outfile, 'w') as out:
        out.write('CDR3 length,fraction cysteines conserved,number of seqs\n')
        for item in cdr_list_all:
            item._calculate_stats()
            out.write(str(item._length) + ',' + str(item._cysteines) + ',' + str(item._number) + '\n')
            data.append([item._length, item._cysteines])
    return data

def plot_conserved_cysteines(infile, outfile, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    data = conserved_cysteines(folders, outfile, regex, frequency)
    arr = np.array(data)
    fig = plt.figure(facecolor = 'white')
    plt.plot(arr[:,0], arr[:,1])
    plt.ylabel('relative frequency of cysteine conservation')
    plt.xlabel('CDR3 length')
    plt.ylim([0,1.2])
    plt.savefig(outfile + '.png')


if __name__=="__main__":
     main(sys.argv) 


                    
                                            
                            
                        
                    
