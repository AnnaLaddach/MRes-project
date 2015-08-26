#Functions to analyse J gene segment distributions of CDR3s with a certain property. 
#A textfile listing IMGT output folders and output file name are taken as input by 'plot_jsegs'. 
#Additionally a fuction to determine whether CDR3s satisfy a certain property must be supplied as an input argument.
#These are included in helper_functions.py.
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID

import csv
import re
import fractions
import numpy as np
import matplotlib.pyplot as plt
import helper_functions


class Cdr:
    def __init__(self):
        self.number = 0 
        self.j_segs = {'IGHJ1':0,'IGHJ2':0,'IGHJ3':0,\
                       'IGHJ4':0, 'IGHJ5':0, 'IGHJ6':0}
        self.fractions = {}
        

    def increment(self, jseg, frequency):
        self.number += frequency
        self.j_segs[jseg] += frequency

    def calculate_stats(self):
        for key in self.j_segs.keys():
            self.fractions[key] = float(self.j_segs[key])/float(self.number)

class Cdr_all_files:
    def __init__(self):
        self.number = 0 
        self.j_segs = {'IGHJ1':0,'IGHJ2':0,'IGHJ3':0,\
                       'IGHJ4':0, 'IGHJ5':0, 'IGHJ6':0}


    def increment(self, jseg, fraction):
        self.number += fraction
        self.j_segs[jseg] += fraction


def analyse_jsegs(folders, outfile, condition, regex = None, frequency = False):
    cdr_all = Cdr_all_files()
    for folder in folders:
        with open(folder[0] + '/5_AA-sequences.txt') as f:
            cdr = Cdr()
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'].find('productive') == 0:
                    if condition(row['CDR3-IMGT']):
                        try:
                            pat = re.compile('IGHJ[^\*]+')
                            info = pat.search(row['J-GENE and allele'])
                            jgene = info.group(0)
                        except:
                            jgene = ''
                        if frequency:
                            pat = re.compile(regex)
                            info = pat.match(row['Sequence ID'])
                            freq = int(info.group(1))
                        else:
                            freq = 1
                        cdr.increment(jgene, freq)
            cdr.calculate_stats()
            for key in cdr.fractions.keys():
                cdr_all.increment(key, cdr.fractions[key] * fractions.Fraction(folder[1]))
    with open(outfile, 'w') as out:
        out.write('IGHJ1,IGHJ2,IGHJ3,IGHJ4,IGHJ5,IGHJ6\n')
        out.write(str(cdr_all.j_segs['IGHJ1']) + ',' + str(cdr_all.j_segs['IGHJ2']) + ',' + str(cdr_all.j_segs['IGHJ3'])\
                              + ','  + str(cdr_all.j_segs['IGHJ4']) + ',' + str(cdr_all.j_segs['IGHJ5']) + ',' + str(cdr_all.j_segs['IGHJ6']) \
                              + '\n')
    jseg_list = [cdr_all.j_segs['IGHJ1'], cdr_all.j_segs['IGHJ2'],cdr_all.j_segs['IGHJ3']\
                              ,cdr_all.j_segs['IGHJ4'],cdr_all.j_segs['IGHJ5'],cdr_all.j_segs['IGHJ6']]
    return jseg_list


def plot_jsegs(infile, outfile, condition, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    data = analyse_jsegs(folders, outfile, condition, regex, frequency)
    ind = np.arange(len(data))
    labels = ['IGHJ1','IGHJ2','IGHJ3','IGHJ4','IGHJ5','IGHJ6']
    fig = plt.figure(facecolor = 'white')
    plt.bar(ind, data, align = 'center')
    plt.ylim([0,0.8])
    plt.ylabel('relative frequency of J-gene')
    plt.xticks(ind, labels, rotation=90, rotation_mode="anchor")
    plt.tick_params(axis = 'x', pad = 30, which = 'both', bottom = 'off', top = 'off')
    plt.tight_layout()
    plt.savefig(outfile + '.png')



                                          
                
                        





