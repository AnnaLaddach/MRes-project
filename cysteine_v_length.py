#Script to analyse and plot the relative frequency of CDR3s containing two or more cysteines in relationship to CDR3 length.
#Naive simulations using the codon based and mammlian protein based frequency of the occurence as well as CDR3 length distributions are also plotted. 
#A textfile listing IMGT output folders and output file name are taken as input. 
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID.

import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import helper_functions
import gen_rand_seqs
import fractions
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.font_manager as fm

def main(argv):
    if len(argv) == 3:
        plot_cysteines_v_length(argv[1], argv[2])
    elif len(argv) == 5:
	plot_cysteines_v_length(argv[1], argv[2], argv[3], argv[4])
    else:
        print 'usage: cysteine_v_length.py IMGT_folder_list outfile '
	print 'optional arguments: frequency_regex frequency(TRUE/FALSE)'
        sys.exit()


class Cdr:
    def __init__(self, length, cysteines, frequency):
        self._length = length
        self._number = 0
        self._cysteines = 0 
        self._fraction = 0
        self._increment(cysteines,frequency)

    def _increment(self,cysteines, frequency):
        self._number += frequency
        self._cysteines += cysteines
        self._fraction = float(self._cysteines)/float(self._number)

class Cdr_all_files:
    def __init__(self, length, fraction_seqs, cysteine_fraction, fraction, number):
        self._length = length
        self._cysteine_fraction = []
        self._fraction_seqs = 0
        self._num_files = 0
        self._number = 0
        self._increment(fraction_seqs, cysteine_fraction, fraction, number)

    def _increment(self, fraction_seqs, cysteine_fraction, fraction, number):
        self._fraction_seqs += fraction_seqs * fraction
        self._cysteine_fraction.append(cysteine_fraction)
        self._number += number

    def _calculate_stats(self):
        self._cysteine_fraction = np.mean(self._cysteine_fraction) 
        


def cysteines_v_length(folders, outfile, regex = None, frequency = False):
    cdr_list_all = []
    for folder in folders:
        with open(folder[0] + '/5_AA-sequences.txt') as f:
            cdr_list = []
            total_seqs = 0
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'] == 'productive':
                    length = len(row['CDR3-IMGT'])
                    if length > 37:
                        print folder
                        print row['Sequence ID']
                    if frequency:
                        pat = re.compile(regex)
                        info = pat.match(row['Sequence ID'])
                        freq = int(info.group(1))
                    else:
                        freq = 1
                    total_seqs += freq
                    cysteines = 0
                    if row['CDR3-IMGT'].count('C') > 1:
                        cysteines = freq
                    notfound = True
                    for item in cdr_list:
                        if length == item._length:
                            item._increment(cysteines, freq)
                            notfound = False
                    if notfound:
                        cdr_list.append(Cdr(length, cysteines, freq))
            for item in cdr_list:
                notfound2 = True
                for item2 in cdr_list_all:
                    if item._length == item2._length:
                        item2._increment(item._number/float(total_seqs),item._fraction, fractions.Fraction(folder[1]), item._number)
                        notfound2 = False
                if notfound2:
                    cdr_list_all.append(Cdr_all_files(item._length, item._number/float(total_seqs),item._fraction, fractions.Fraction(folder[1]), item._number))
    cdr_list_all.sort(key = lambda x: x._length)
    data = []
    with open(outfile, 'w') as out:
        out.write('CDR3 length,fraction of sequences,'\
                   + 'fraction of CDR3s with 2 or more cysteines,number of sequences\n')
        for item in cdr_list_all:
            item._calculate_stats()
            out.write(str(item._length) + ',' + str(item._fraction_seqs) + ',' + str(item._cysteine_fraction) + ',' + str(item._number) + '\n')
            data.append([item._length,item._fraction_seqs,item._cysteine_fraction])
    return data

def plot_cysteines_v_length(infile, outfile, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    font = fm.FontProperties(family = 'Sans-Serif', size = 20)
    font2 = fm.FontProperties(family = 'Sans-Serif', size = 25)
    folders = helper_functions.create_list(infile)
    data = cysteines_v_length(folders, outfile, regex, frequency)
    values = range(9)
    data1 = gen_rand_seqs.generate_random_seqs_2(range(1,40),1000000)
    data2 = gen_rand_seqs.generate_random_seqs_3(range(1,40),1000000)
    jet = cm = plt.get_cmap('CMRmap') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    colorVal = scalarMap.to_rgba(values[4])
    colorVal1 = scalarMap.to_rgba(values[0])
    colorVal2 = scalarMap.to_rgba(values[2])
    arr = np.array(data)
    arr1 = np.array(data1)
    arr2 = np.array(data2)
    fig, ax1 = plt.subplots()
    fig.set_size_inches(10, 8)
    ax1.plot(arr[:,0], arr[:,1], ls = '--', color = colorVal2)
    ax1.scatter(arr[:,0], arr[:,1], color = colorVal2)
    plt.ylim(0,0.14)
    ax1.set_xlabel('CDR3 length (amino acids)')
    ax1.set_ylabel('relative frequency (n > 100)', color = colorVal2)
    for t1 in ax1.get_yticklabels():
        t1.set_color(colorVal2)
        t1.set_fontproperties(font)
    ax2 = ax1.twinx()
    ax2.plot(arr[:,0],arr[:,2], color = colorVal)
    ax2.scatter(arr[:,0],arr[:,2], color = colorVal)
    ax2.plot(arr1[:,0],arr1[:,1], color = colorVal, linestyle = ':', label = 'simulation 1')
    ax2.plot(arr2[:,0],arr2[:,1], color = colorVal, linestyle = '-.', label = 'simulation 2' )
    plt.ylim(0,0.5)
    plt.xlim(0,40)
    plt.legend(prop = font2)
    ax2.set_ylabel('relative frequency of CDR3s containing two or more cysteines', color = colorVal)
    for t2 in ax2.get_yticklabels():
        t2.set_color(colorVal)
        t2.set_fontproperties(font)
    for t in ax1.get_xticklabels():
        t.set_fontproperties(font)
    for t in ax2.get_xticklabels():
        t.set_fontproperties(font)
    plt.savefig(outfile + '.png')

if __name__=="__main__":
     main(sys.argv) 
 
    

        
