#Functions to analyse and plot 5' and 3' N-nucleotide use as a function of CDR3 length in CDRs with a defined property.
#A textfile listing IMGT output folders and output file name are taken as input for 'plot_n_nucleotides'. 
#Additionally a fuction to determine whether CDR3s satisfy a certain property must be supplied as an input argument.
#These are included in helper_functions.py.
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID


import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import helper_functions
import fractions

class Cdr:
    def __init__(self, length, n1, n2, frequency):
        self._length = length
        self._number = 0
        self._n1 = 0
        self._n2 = 0
        self._mean_n1 = 0
        self._mean_n2 = 0
        self._increment(n1, n2, frequency)

    def _increment(self, n1, n2, frequency):
        self._n1 += n1 * frequency
        self._n2 += n2 * frequency
        self._number += frequency

    def _calculate_stats(self):
        self._mean_n1 = self._n1 / self._number
        self._mean_n2 = self._n2 / self._number

class Cdr_all:
    def __init__(self, length, n1, n2):
        self._length = length
        self._n1 = []
        self._n2 = []
        self._number = 0
        self._increment(n1, n2)

    def _increment(self, n1, n2):
        self._n1.append(n1)
        self._n2.append(n2)
        self._number += 1

    def _calculate_stats(self):
        self._n1 = np.mean(self._n1)
        self._n2 = np.mean(self._n2) 
    

def n_nucleotides(folders, outfile, condition, regex = None, frequency = False):
    cdr_list_all = []
    for folder in folders:
        with open(folder[0] + '/6_Junction.txt') as f:
            cdr_list = []
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'].find('productive') == 0:
                    if condition(row['CDR3-IMGT (AA)']):
                        if row['N-REGION-nt nb'] == '0':
                            n1 = int(row['N1-REGION-nt nb'])
                            n2 = int(row['N2-REGION-nt nb'])
                            length = len(row['CDR3-IMGT (AA)'])
                            if frequency:
                                pat = re.compile(regex)
                                info = pat.match(row['Sequence ID'])
                                freq = int(info.group(1))
                            else:
                                freq = 1
                            notfound = True
                            for item in cdr_list:
                                if length == item._length:
                                    item._increment(n1, n2, freq)
                                    notfound = False
                            if notfound:
                                cdr_list.append(Cdr(length, n1, n2, freq))
            for item in cdr_list:
                item._calculate_stats()
                not_found2 = True
                for item2 in cdr_list_all:
                    if item._length == item2._length:
                        item2._increment(item._mean_n1, item._mean_n2)
                        not_found2 = False
                if not_found2:
                    cdr_list_all.append(Cdr_all(item._length, item._mean_n1, item._mean_n2))
    cdr_list_all.sort(key = lambda x: x._length)
    data = []
    with open(outfile, 'w') as out:
        out.write('CDR3 length,mean n1 nucleotides,mean n2 nucleotides\n')
        for item in cdr_list_all:
            item._calculate_stats()
            out.write(str(item._length) + ',' + str(item._n1) + ',' + str(item._n2) +'\n')
            data.append([item._length, item._n1, item._n2])
    return data


def plot_n_nucleotides(infile, outfile, condition, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    data = n_nucleotides(folders, outfile, condition, regex, frequency)
    arr = np.array(data)
    ind = np.arange(len(data))
    width = 0.35
    fig, ax = plt.subplots()
    n1_bars = ax.bar(ind, arr[:,1], width, color = 'b')
    n2_bars = ax.bar(ind + width, arr[:,2], width, color = 'r')
    ax.set_ylabel('number of n nucleotides')
    ax.set_ylim([0,60])
    ax.set_xlabel('CDR3 length')
    ax.set_xticks(ind + width)
    ax.set_xticklabels([int(i) for i in arr[:,0]])
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=90 )
    ax.legend((n1_bars[0], n2_bars[0]),('n1', 'n2'))
    plt.savefig(outfile + '.png')                                                                                 
                
