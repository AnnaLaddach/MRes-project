#Functions to analyse V gene segment distributions of CDR3s with a certain property. 
#A textfile listing IMGT output folders and output file name are taken as input by 'plot_vsegs'. 
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
        self.v_segs = {'IGHV1':0,'IGHV2':0,'IGHV3':0,'IGHV4':0,'IGHV5':0,\
                       'IGHV6':0,'IGHV7':0}
        self.fractions = {}
        

    def increment(self, vseg, frequency):
        self.number += frequency
        self.v_segs[vseg] += frequency

    def calculate_stats(self):
        for key in self.v_segs.keys():
            self.fractions[key] = float(self.v_segs[key])/float(self.number)

class Cdr_all_files:
    def __init__(self):
        self.number = 0 
        self.v_segs = {'IGHV1':[],'IGHV2':[],'IGHV3':[],'IGHV4':[],'IGHV5':[],\
                       'IGHV6':[],'IGHV7':[]}
        self.v_segs_mean = {'IGHV1':0,'IGHV2':0,'IGHV3':0,'IGHV4':0,'IGHV5':0,\
                       'IGHV6':0,'IGHV7':0}
        self.v_segs_sd = {'IGHV1':0,'IGHV2':0,'IGHV3':0,'IGHV4':0,'IGHV5':0,\
                       'IGHV6':0,'IGHV7':0}

    def increment(self, vseg, fraction):
        self.v_segs[vseg].append(fraction)

    def increment_number(self):
        self.number += 1
        
    def calculate_stats(self):
        for key in self.v_segs.keys():
            while len(self.v_segs[key]) < self.number:
                self.v_segs[key].append(0)
            self.v_segs_mean[key] = np.mean(self.v_segs[key])
            self.v_segs_sd[key] = np.std(self.v_segs[key])
            


def analyse_vsegs(folders, outfile, condition, regex = None, frequency = False):
    cdr_all = Cdr_all_files()
    for folder in folders:
        with open(folder[0] + '/5_AA-sequences.txt') as f:
            cdr = Cdr()
            cdr_all.increment_number()
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'].find('productive') == 0:
                    if condition(row['CDR3-IMGT']):
                        vgene = helper_functions.find_v_gene(row['V-GENE and allele'])
                        if vgene:
                            if frequency:
                                pat = re.compile(regex)
                                info = pat.match(row['Sequence ID'])
                                freq = int(info.group(1))
                            else:
                                freq = 1
                            cdr.increment(vgene, freq)
            cdr.calculate_stats()
            for key in cdr.fractions.keys():
                cdr_all.increment(key, cdr.fractions[key])
    print cdr_all.number
    cdr_all.calculate_stats()
    with open(outfile, 'w') as out:
        out.write('IGHV1,IGHV2,IGHV3,IGHV4,IGHV5,IGHV6,IGHV7\n')
        out.write(str(cdr_all.v_segs_mean['IGHV1']) + ',' + str(cdr_all.v_segs_mean['IGHV2']) + ',' + str(cdr_all.v_segs_mean['IGHV3']) + ',' + str(cdr_all.v_segs_mean['IGHV4'])\
                              + ',' + str(cdr_all.v_segs_mean['IGHV5']) + ',' + str(cdr_all.v_segs_mean['IGHV6']) + ',' + str(cdr_all.v_segs_mean['IGHV7'])\
                              +',' + str(cdr_all.v_segs_sd['IGHV1']) + ',' + str(cdr_all.v_segs_sd['IGHV2']) + ',' + str(cdr_all.v_segs_sd['IGHV3']) + ',' + str(cdr_all.v_segs_sd['IGHV4'])\
                              + ',' + str(cdr_all.v_segs_sd['IGHV5']) + ',' + str(cdr_all.v_segs_sd['IGHV6']) + ',' + str(cdr_all.v_segs_sd['IGHV7']) + '\n')
    vseg_list = [[cdr_all.v_segs_mean['IGHV1'],cdr_all.v_segs_mean['IGHV2'],cdr_all.v_segs_mean['IGHV3'],cdr_all.v_segs_mean['IGHV4']\
                              ,cdr_all.v_segs_mean['IGHV5'],cdr_all.v_segs_mean['IGHV6'],cdr_all.v_segs_mean['IGHV7']],\
                 [cdr_all.v_segs_sd['IGHV1'],cdr_all.v_segs_sd['IGHV2'],cdr_all.v_segs_sd['IGHV3'],cdr_all.v_segs_sd['IGHV4']\
                              ,cdr_all.v_segs_sd['IGHV5'],cdr_all.v_segs_sd['IGHV6'],cdr_all.v_segs_sd['IGHV7']]]
    return vseg_list



def plot_vsegs(infile, outfile, condition, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    data = analyse_vsegs(folders, outfile, condition, regex, frequency)
    ind = np.arange(len(data[0]))
    labels = ['IGHV1','IGHV2','IGHV3','IGHV4','IGHV5','IGHV6','IGHV7']
    fig = plt.figure(facecolor = 'white')
    plt.bar(ind, data[0], yerr = data[1], align = 'center')
    plt.ylim([0,0.8])
    plt.ylabel('relative frequency of V-gene')
    plt.xticks(ind, labels, rotation=90, rotation_mode="anchor")
    plt.tick_params(axis = 'x', pad = 30, which = 'both', bottom = 'off', top = 'off')
    plt.tight_layout()
    plt.savefig(outfile + '.png')                  
                
                        





