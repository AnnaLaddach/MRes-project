#Script to analyse and plot the positioning of CDR3 5' and 3' cysteines in relationtionship to CDR3 length. 
#A textfile listing IMGT output folders and output file name are taken as input. 
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID.

import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import helper_functions
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.font_manager as fm

def main(argv):
    if len(argv) == 3:
        plot_cysteine_positions(argv[1], argv[2])
    elif len(argv) == 5:
	plot_cysteine_positions(argv[1], argv[2], argv[3], argv[4])
    else:
        print 'usage: positioning_of_cysteines.py IMGT_folder_list outfile '
	print 'optional arguments: frequency_regex frequency(TRUE/FALSE)'
        sys.exit()

class Cdr:
    def __init__(self, length, start, end, frequency):
        self._length = length
        self._number = 0
        self._start = []
        self._start_mean = 0
        self._start_stdev = 0
        self._end = []
        self._end_mean =0
        self._end_stdev = 0
        self._increment(start, end, frequency)

    def _increment(self, start, end, frequency):
        self._number += frequency
        for i in range(frequency):
            self._start.append(start)
            self._end.append(end)
            
    def _calculate_stats(self):
        self._start_mean = np.mean(self._start)
        self._end_mean = np.mean(self._end)
        self._start_stdev = np.std(self._start)
        self._end_stdev = np.std(self._end)

class Cdr_all_files:
    def __init__(self, length, start_mean, end_mean, start_stdev, end_stdev, number):
        self._length = length
        self._number = 0
        self._num_files = 0
        self._start_mean = 0
        self._start_stdev = 0
        self._end_mean = 0
        self._end_stdev = 0
        self._increment(start_mean, end_mean, start_stdev, end_stdev, number)

    def _increment(self, start_mean, end_mean, start_stdev, end_stdev,number):
        self._num_files += 1
        self._start_mean += start_mean
        self._start_stdev += start_stdev
        self._end_mean += end_mean
        self._end_stdev += end_stdev
        self._number += number

    def _calculate_stats(self):
        self._start_mean = self._start_mean/float(self._num_files)
        self._end_mean = self._end_mean/float(self._num_files)
        self._start_stdev = self._start_stdev/float(self._num_files)
        self._end_stdev = self._end_stdev/float(self._num_files)
          
    
def find_cysteine_positioning(folders, outfile, regex = None, frequency = False):
    cdr_list_all = []
    for folder in folders:
        with open(folder[0] + '/5_AA-sequences.txt') as f:  
            cdr_list = []
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'].find('productive') == 0:
                    if row['CDR3-IMGT'].count('C') > 1:
                        if frequency:
                            pat = re.compile(regex)
                            info = pat.match(row['Sequence ID'])
                            freq = int(info.group(1))
                        else:
                            freq = 1
                        length = len(row['CDR3-IMGT'])
                        notfound = True
                        for item in cdr_list:
                            if length == item._length:
                                item._increment((row['CDR3-IMGT'].find('C') + 1)/float(length), (length - row['CDR3-IMGT'].rfind('C'))/float(length), freq)
                                notfound = False
                        if notfound:
                            cdr_list.append(Cdr(length, (row['CDR3-IMGT'].find('C') + 1)/float(length), (length - row['CDR3-IMGT'].rfind('C'))/float(length), freq))
            for item in cdr_list:
                item._calculate_stats()
                notfound2 = True
                for item2 in cdr_list_all:
                    if item._length == item2._length:
                        item2._increment(item._start_mean, item._end_mean, item._start_stdev, item._end_stdev, item._number)
                        notfound2 = False
                if notfound2:
                    cdr_list_all.append(Cdr_all_files(item._length, item._start_mean, item._end_mean, item._start_stdev, item._end_stdev, item._number))         
    cdr_list_all.sort(key = lambda x: x._length)
    data = []
    with open(outfile, 'w') as out:
        out.write('CDR3 length,mean start position,stdev start position,mean stop position,stdev stop position,n sequences\n')
        for item in cdr_list_all:
            item._calculate_stats()
            out.write(str(item._length) + ',' + str(item._start_mean) + ',' + str(item._start_stdev)\
                      + ',' + str(item._end_mean) + ',' + str(item._end_stdev) + ',' + str(item._number) + '\n')
            data.append([item._length, item._start_mean, item._start_stdev, item._end_mean, item._end_stdev])
    return data


def plot_cysteine_positions(infile, outfile, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    data = find_cysteine_positioning(folders, outfile, regex, frequency)
    arr = np.array(data)
    n = len(data)
    ind = np.arange(n)
    font = fm.FontProperties(family = 'Sans-Serif', size = 20)
    font2 = fm.FontProperties(family = 'Sans-Serif', size = 25)
    font1 = fm.FontProperties(family = 'Sans-Serif', size = 15)
    values = range(9)
    c = cm = plt.get_cmap('CMRmap') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=c)
    width = 0.35
    fig, ax = plt.subplots(figsize = (20,10))
    start_bars = ax.bar(ind + 0.15, arr[:,1], width, color =  scalarMap.to_rgba(values[2]))
    stop_bars = ax.bar(ind + 0.5, arr[:,3], width, color= scalarMap.to_rgba(values[4]))
    ax.set_ylabel('mean distance of cysteine from start/end of CDR3 as a fraction of CDR3 length', fontproperties = font1)
    ax.set_xlabel('CDR3 length (amino acids)', fontproperties = font2)
    ax.set_xticks(ind + 0.5)
    ax.set_xticklabels([int(i) for i in arr[:,0]])
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontproperties(font)
    ax.legend((start_bars[0], stop_bars[0]),('5\' cysteine', '3\' cysteine'), prop = font2)
    plt.savefig(outfile + '.png')

if __name__=="__main__":
     main(sys.argv) 

