#Functions to analyse D gene segment distributions of CDR3s with a certain property. 
#A textfile listing IMGT output folders and output file name are taken as input by 'plot_dsegs'. 
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
        self.d_segs = {'IGHD1-1':0,'IGHD1-7':0,'IGHD1-14':0,'IGHD1-20':0,'IGHD1-26':0,\
                       'IGHD2-2':0,'IGHD2-8':0, 'IGHD2-15':0, 'IGHD2-21':0,\
                       'IGHD3-3':0,'IGHD3-9':0,'IGHD3-10':0,'IGHD3-16':0,'IGHD3-22':0,
                       'IGHD4-4':0,'IGHD4-11':0,'IGHD4-17':0,'IGHD4-23':0,\
                       'IGHD5-5':0,'IGHD5-12':0,'IGHD5-18':0,'IGHD5-24':0,\
                       'IGHD6-6':0,'IGHD6-13':0,'IGHD6-19':0,'IGHD6-25':0,'IGHD7-27':0, '': 0}
        self.d_segs_nc = {'IGHD1-1':0,'IGHD1-7':0,'IGHD1-14':0,'IGHD1-20':0,'IGHD1-26':0,\
                       'IGHD2-2':0,'IGHD2-8':0, 'IGHD2-15':0, 'IGHD2-21':0,\
                       'IGHD3-3':0,'IGHD3-9':0,'IGHD3-10':0,'IGHD3-16':0,'IGHD3-22':0,
                       'IGHD4-4':0,'IGHD4-11':0,'IGHD4-17':0,'IGHD4-23':0,\
                       'IGHD5-5':0,'IGHD5-12':0,'IGHD5-18':0,'IGHD5-24':0,\
                       'IGHD6-6':0,'IGHD6-13':0,'IGHD6-19':0,'IGHD6-25':0,'IGHD7-27':0, '': 0}
        self.fractions = {}
        self.fractions_nc = {}
        
    def increment_nc(self, dseg, frequency):
        self.number += frequency
        self.d_segs_nc[dseg] += frequency

    def increment(self, dseg, frequency):
        self.number += frequency
        self.d_segs[dseg] += frequency

    def calculate_stats(self):
        for key in self.d_segs.keys():
            self.fractions[key] = float(self.d_segs[key])/float(self.number)
        for key in self.d_segs_nc.keys():
            self.fractions_nc[key] = float(self.d_segs_nc[key])/float(self.number)


class Cdr_all_files:
    def __init__(self):
        self.d_segs = {'IGHD1-1':0,'IGHD1-7':0,'IGHD1-14':0,'IGHD1-20':0,'IGHD1-26':0,\
                       'IGHD2-2':0,'IGHD2-8':0, 'IGHD2-15':0, 'IGHD2-21':0,\
                       'IGHD3-3':0,'IGHD3-9':0,'IGHD3-10':0,'IGHD3-16':0,'IGHD3-22':0,
                       'IGHD4-4':0,'IGHD4-11':0,'IGHD4-17':0,'IGHD4-23':0,\
                       'IGHD5-5':0,'IGHD5-12':0,'IGHD5-18':0,'IGHD5-24':0,\
                       'IGHD6-6':0,'IGHD6-13':0,'IGHD6-19':0,'IGHD6-25':0,'IGHD7-27':0, '': 0}
        self.d_segs_nc = {'IGHD1-1':0,'IGHD1-7':0,'IGHD1-14':0,'IGHD1-20':0,'IGHD1-26':0,\
                       'IGHD2-2':0,'IGHD2-8':0, 'IGHD2-15':0, 'IGHD2-21':0,\
                       'IGHD3-3':0,'IGHD3-9':0,'IGHD3-10':0,'IGHD3-16':0,'IGHD3-22':0,
                       'IGHD4-4':0,'IGHD4-11':0,'IGHD4-17':0,'IGHD4-23':0,\
                       'IGHD5-5':0,'IGHD5-12':0,'IGHD5-18':0,'IGHD5-24':0,\
                       'IGHD6-6':0,'IGHD6-13':0,'IGHD6-19':0,'IGHD6-25':0,'IGHD7-27':0, '': 0}

    def increment(self, dseg, fraction):
        self.d_segs[dseg] += fraction

    def increment_nc(self, dseg, fraction):
        self.d_segs_nc[dseg] += fraction

    def alleles(self):
        self.d_segs['IGHD4-4/4-11'] = self.d_segs['IGHD4-4'] + self.d_segs['IGHD4-11']
        self.d_segs['IGHD5-5/5-18'] = self.d_segs['IGHD5-5'] + self.d_segs['IGHD5-18']
        self.d_segs_nc['IGHD4-4/4-11'] = self.d_segs_nc['IGHD4-4'] + self.d_segs_nc['IGHD4-11']
        self.d_segs_nc['IGHD5-5/5-18'] = self.d_segs_nc['IGHD5-5'] + self.d_segs_nc['IGHD5-18']



def analyse_dsegs(folders, outfile, regex = None, frequency = False):
    cdr_all = Cdr_all_files()
    for folder in folders:
        with open(folder[0] + '/5_AA-sequences.txt') as f:
            cdr = Cdr()
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'].find('productive') == 0:
                    if row['CDR3-IMGT'].count('C') > 1:
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
                        if helper_functions.non_canonical(row['CDR3-IMGT']):
                            cdr.increment_nc(dgene, freq)
                        else:
                            cdr.increment(dgene, freq) 
            cdr.calculate_stats()
            for key in cdr.fractions.keys():
                cdr_all.increment(key, cdr.fractions[key] * fractions.Fraction(folder[1]))
                cdr_all.increment_nc(key, cdr.fractions_nc[key] * fractions.Fraction(folder[1]))  
    cdr_all.alleles()
    with open(outfile, 'w') as out:
        out.write('IGHD1-1nc,IGHD1-7nc,IGHD1-14nc,IGHD1-20nc,IGHD1-26nc,\
            IGHD2-2nc,IGHD2-8nc,IGHD2-15nc,IGHD2-21nc,\
            IGHD3-3nc,IGHD3-9nc,IGHD3-10nc,IGHD3-16nc,IGHD3-22nc,\
            IGHD4-4nc,IGHD4-11nc,IGHD4-17nc,IGHD4-23nc,\
            IGHD5-5nc,IGHD5-12nc,IGHD5-18nc,IGHD5-24nc,\
            IGHD6-6nc,IGHD6-13nc,IGHD6-19nc,IGHD6-25nc,IGHD7-27nc,nonenc\
            IGHD1-1,IGHD1-7,IGHD1-14,IGHD1-20,IGHD1-26,\
            IGHD2-2,IGHD2-8,IGHD2-15,IGHD2-21,\
            IGHD3-3,IGHD3-9,IGHD3-10,IGHD3-16,IGHD3-22,\
            IGHD4-4,IGHD4-11,IGHD4-17,IGHD4-23,\
            IGHD5-5,IGHD5-12,IGHD5-18,IGHD5-24,\
            IGHD6-6,IGHD6-13,IGHD6-19,IGHD6-25,IGHD7-27,none\n')
        out.write(str(cdr_all.d_segs_nc['IGHD1-1']) + ',' + str(cdr_all.d_segs_nc['IGHD1-7']) + ',' + str(cdr_all.d_segs_nc['IGHD1-14']) + ',' + str(cdr_all.d_segs_nc['IGHD1-20']) + ',' + str(cdr_all.d_segs_nc['IGHD1-26'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD2-2']) + ',' + str(cdr_all.d_segs_nc['IGHD2-8']) + ',' + str(cdr_all.d_segs_nc['IGHD2-15']) + ',' + str(cdr_all.d_segs_nc['IGHD2-21']) \
                              + ',' + str(cdr_all.d_segs_nc['IGHD3-3']) + ',' + str(cdr_all.d_segs_nc['IGHD3-9'])  + ',' + str(cdr_all.d_segs_nc['IGHD3-10'])  + ',' + str(cdr_all.d_segs_nc['IGHD3-16'])  + ',' + str(cdr_all.d_segs_nc['IGHD3-22'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD4-4/4-11']) + ',' + str(cdr_all.d_segs_nc['IGHD4-17']) + ',' + str(cdr_all.d_segs_nc['IGHD4-23'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD5-5/5-18']) + ',' + str(cdr_all.d_segs_nc['IGHD5-12']) + ',' + str(cdr_all.d_segs_nc['IGHD5-24'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD6-6']) + ',' + str(cdr_all.d_segs_nc['IGHD6-13']) + ',' + str(cdr_all.d_segs_nc['IGHD6-19']) + ',' + str(cdr_all.d_segs_nc['IGHD6-25'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD7-27']) + ',' + str(cdr_all.d_segs_nc['']) + ','\
                              + str(cdr_all.d_segs['IGHD1-1']) + ',' + str(cdr_all.d_segs['IGHD1-7']) + ',' + str(cdr_all.d_segs['IGHD1-14']) + ',' + str(cdr_all.d_segs['IGHD1-20']) + ',' + str(cdr_all.d_segs['IGHD1-26'])\
                              + ',' + str(cdr_all.d_segs['IGHD2-2']) + ',' + str(cdr_all.d_segs['IGHD2-8']) + ',' + str(cdr_all.d_segs['IGHD2-15']) + ',' + str(cdr_all.d_segs['IGHD2-21']) \
                              + ',' + str(cdr_all.d_segs['IGHD3-3']) + ',' + str(cdr_all.d_segs['IGHD3-9'])  + ',' + str(cdr_all.d_segs['IGHD3-10'])  + ',' + str(cdr_all.d_segs['IGHD3-16'])  + ',' + str(cdr_all.d_segs['IGHD3-22'])\
                              + ',' + str(cdr_all.d_segs['IGHD4-4/4-11']) + ',' + str(cdr_all.d_segs['IGHD4-17']) + ',' + str(cdr_all.d_segs['IGHD4-23'])\
                              + ',' + str(cdr_all.d_segs['IGHD5-5/5-18']) + ',' + str(cdr_all.d_segs['IGHD5-12']) + ',' + str(cdr_all.d_segs['IGHD5-24'])\
                              + ',' + str(cdr_all.d_segs['IGHD6-6']) + ',' + str(cdr_all.d_segs['IGHD6-13']) + ',' + str(cdr_all.d_segs['IGHD6-19']) + ',' + str(cdr_all.d_segs['IGHD6-25'])\
                              + ',' + str(cdr_all.d_segs['IGHD7-27']) + ',' + str(cdr_all.d_segs['']) + '\n')
    dseg_list_nc = [cdr_all.d_segs_nc['IGHD1-1'],cdr_all.d_segs_nc['IGHD1-7'],cdr_all.d_segs_nc['IGHD1-14'],cdr_all.d_segs_nc['IGHD1-20'],cdr_all.d_segs_nc['IGHD1-26']\
                              ,cdr_all.d_segs_nc['IGHD2-2'],cdr_all.d_segs_nc['IGHD2-8'],cdr_all.d_segs_nc['IGHD2-15'],cdr_all.d_segs_nc['IGHD2-21'] \
                              ,cdr_all.d_segs_nc['IGHD3-3'],cdr_all.d_segs_nc['IGHD3-9'],cdr_all.d_segs_nc['IGHD3-10'],cdr_all.d_segs_nc['IGHD3-16'],cdr_all.d_segs_nc['IGHD3-22']\
                              ,cdr_all.d_segs_nc['IGHD4-4/4-11'],cdr_all.d_segs_nc['IGHD4-17'],cdr_all.d_segs_nc['IGHD4-23']\
                              ,cdr_all.d_segs_nc['IGHD5-5/5-18'],cdr_all.d_segs_nc['IGHD5-12'], cdr_all.d_segs_nc['IGHD5-24']\
                              ,cdr_all.d_segs_nc['IGHD6-6'],cdr_all.d_segs_nc['IGHD6-13'],cdr_all.d_segs_nc['IGHD6-19'],cdr_all.d_segs_nc['IGHD6-25'],cdr_all.d_segs_nc['IGHD7-27'],cdr_all.d_segs_nc['']]
    dseg_list = [cdr_all.d_segs['IGHD1-1'],cdr_all.d_segs['IGHD1-7'],cdr_all.d_segs['IGHD1-14'],cdr_all.d_segs['IGHD1-20'],cdr_all.d_segs['IGHD1-26']\
                              ,cdr_all.d_segs['IGHD2-2'],cdr_all.d_segs['IGHD2-8'],cdr_all.d_segs['IGHD2-15'],cdr_all.d_segs['IGHD2-21'] \
                              ,cdr_all.d_segs['IGHD3-3'],cdr_all.d_segs['IGHD3-9'],cdr_all.d_segs['IGHD3-10'],cdr_all.d_segs['IGHD3-16'],cdr_all.d_segs['IGHD3-22']\
                              ,cdr_all.d_segs['IGHD4-4/4-11'],cdr_all.d_segs['IGHD4-17'],cdr_all.d_segs['IGHD4-23']\
                              ,cdr_all.d_segs['IGHD5-5/5-18'],cdr_all.d_segs['IGHD5-12'], cdr_all.d_segs['IGHD5-24']\
                              ,cdr_all.d_segs['IGHD6-6'],cdr_all.d_segs['IGHD6-13'],cdr_all.d_segs['IGHD6-19'],cdr_all.d_segs['IGHD6-25'],cdr_all.d_segs['IGHD7-27'],cdr_all.d_segs['']]
    return [dseg_list_nc, dseg_list]



def plot_dsegs(infile, outfile, condition, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    data = analyse_dsegs(folders, outfile, condition, regex, frequency)
    print len(data)
    ind = np.arange(len(data))
    labels = ['IGHD1-1','IGHD1-7','IGHD1-14','IGHD1-20','IGHD1-26',\
            'IGHD2-2','IGHD2-8','IGHD2-15','IGHD2-21',\
            'IGHD3-3','IGHD3-9','IGHD3-10','IGHD3-16','IGHD3-22',\
            'IGHD4-4/4-11','IGHD4-17','IGHD4-23',\
            'IGHD5-5/5-18','IGHD5-12','IGHD5-24',\
            'IGHD6-6','IGHD6-13','IGHD6-19','IGHD6-25','IGHD7-27','none']
    print len(ind)
    fig = plt.figure(facecolor = 'white')
    plt.bar(ind, data)
    plt.ylim([0,0.5])
    plt.ylabel('relative frequency of D-gene')
    plt.xticks(ind, labels, rotation=90, rotation_mode="anchor")
    plt.tick_params(axis = 'x', pad = 40, which = 'both', bottom = 'off', top = 'off')
    plt.tight_layout()
    plt.savefig(outfile + '.png')


                                          
                
                        





