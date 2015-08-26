#Functions to analyse macaque D gene segment distributions of CDR3s with a certain property. 
#A textfile listing IMGT output folders and output file name are taken as input by 'plot_dsegs'. 
#Additionally a fuction to determine whether CDR3s satisfy a certain property must be supplied as an input argument.
#These are included in helper_functions.py.
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID


import csv
import re
import math
import fractions
import numpy as np
import matplotlib.pyplot as plt
import helper_functions
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.font_manager as fm


class Cdr:
    def __init__(self):
        self.number = 0 
        self.d_segs = {'IGHD1-1':0,'IGHD1-2':0,'IGHD1-3':0,'IGHD1-4':0,'IGHD1-5':0, 'IGHD1-6':0, 'IGHD1-7':0, 'IGHD1-8':0,\
                       'IGHD2-1':0,'IGHD2-2':0, 'IGHD2-3':0, 'IGHD2-4':0, 'IGHD2-5':0, 'IGHD2-6':0,\
                       'IGHD3-1':0,'IGHD3-2':0,'IGHD3-3':0,'IGHD3-4':0,'IGHD4-1':0,
                       'IGHD4-2':0,'IGHD4-3':0,'IGHD4-4':0,'IGHD5-1':0,\
                       'IGHD5-2':0,'IGHD5-3':0,'IGHD6-1':0,'IGHD6-2':0,\
                       'IGHD6-3':0,'IGHD6-4':0,'IGHD6-5':0,'IGHD6-6':0,'IGHD7-1':0, '': 0}
        self.d_segs_nc = {'IGHD1-1':0,'IGHD1-2':0,'IGHD1-3':0,'IGHD1-4':0,'IGHD1-5':0, 'IGHD1-6':0, 'IGHD1-7':0, 'IGHD1-8':0,\
                       'IGHD2-1':0,'IGHD2-2':0, 'IGHD2-3':0, 'IGHD2-4':0, 'IGHD2-5':0, 'IGHD2-6':0,\
                       'IGHD3-1':0,'IGHD3-2':0,'IGHD3-3':0,'IGHD3-4':0,'IGHD4-1':0,
                       'IGHD4-2':0,'IGHD4-3':0,'IGHD4-4':0,'IGHD5-1':0,\
                       'IGHD5-2':0,'IGHD5-3':0,'IGHD6-1':0,'IGHD6-2':0,\
                       'IGHD6-3':0,'IGHD6-4':0,'IGHD6-5':0,'IGHD6-6':0,'IGHD7-1':0, '': 0}
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
        self.d_segs = {'IGHD1-1':0,'IGHD1-2':0,'IGHD1-3':0,'IGHD1-4':0,'IGHD1-5':0, 'IGHD1-6':0, 'IGHD1-7':0, 'IGHD1-8':0,\
                       'IGHD2-1':0,'IGHD2-2':0, 'IGHD2-3':0, 'IGHD2-4':0, 'IGHD2-5':0, 'IGHD2-6':0,\
                       'IGHD3-1':0,'IGHD3-2':0,'IGHD3-3':0,'IGHD3-4':0,'IGHD4-1':0,
                       'IGHD4-2':0,'IGHD4-3':0,'IGHD4-4':0,'IGHD5-1':0,\
                       'IGHD5-2':0,'IGHD5-3':0,'IGHD6-1':0,'IGHD6-2':0,\
                       'IGHD6-3':0,'IGHD6-4':0,'IGHD6-5':0,'IGHD6-6':0,'IGHD7-1':0, '': 0}
        self.d_segs_nc = {'IGHD1-1':0,'IGHD1-2':0,'IGHD1-3':0,'IGHD1-4':0,'IGHD1-5':0, 'IGHD1-6':0, 'IGHD1-7':0, 'IGHD1-8':0,\
                       'IGHD2-1':0,'IGHD2-2':0, 'IGHD2-3':0, 'IGHD2-4':0, 'IGHD2-5':0, 'IGHD2-6':0,\
                       'IGHD3-1':0,'IGHD3-2':0,'IGHD3-3':0,'IGHD3-4':0,'IGHD4-1':0,
                       'IGHD4-2':0,'IGHD4-3':0,'IGHD4-4':0,'IGHD5-1':0,\
                       'IGHD5-2':0,'IGHD5-3':0,'IGHD6-1':0,'IGHD6-2':0,\
                       'IGHD6-3':0,'IGHD6-4':0,'IGHD6-5':0,'IGHD6-6':0,'IGHD7-1':0, '': 0}

    def increment(self, dseg, fraction):
        self.d_segs[dseg] += fraction

    def increment_nc(self, dseg, fraction):
        self.d_segs_nc[dseg] += fraction


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
                        if helper_functions.spacing_4(row['CDR3-IMGT']):
                            cdr.increment(dgene, freq)
                        else:
                            cdr.increment_nc(dgene, freq)
                            
            cdr.calculate_stats()
            for key in cdr.fractions.keys():
                cdr_all.increment(key, cdr.fractions[key] * fractions.Fraction(folder[1]))
                cdr_all.increment_nc(key, cdr.fractions_nc[key] * fractions.Fraction(folder[1]))
                
    #cdr_all.alleles()
    with open(outfile, 'w') as out:
        out.write('IGHD1-1nc,IGHD1-2nc,IGHD1-3nc,IGHD1-4nc,IGHD1-5nc,IGHD1-6nc,IGHD1-7nc,IGHD1-8nc,\
            IGHD2-1nc,IGHD2-2nc,IGHD2-3nc,IGHD2-4nc,IGHD2-5nc,IGHD2-6nc,\
            IGHD3-1nc,IGHD3-2nc,IGHD3-3nc,IGHD3-4nc,IGHD4-1nc,\
            IGHD4-2nc,IGHD4-3nc,IGHD4-4nc,IGHD5-1nc,\
            IGHD5-2nc,IGHD5-3nc,IGHD6-1nc,IGHD6-2nc,\
            IGHD6-3nc,IGHD6-4nc,IGHD6-5nc,IGHD6-6nc,IGHD7-1nc,nonenc,\
            IGHD1-1,IGHD1-2,IGHD1-3,IGHD1-4,IGHD1-5,IGHD1-6,IGHD1-7,IGHD1-8,\
            IGHD2-1,IGHD2-2,IGHD2-3,IGHD2-4,IGHD2-5,IGHD2-6,\
            IGHD3-1,IGHD3-2,IGHD3-3,IGHD3-4,IGHD4-1,\
            IGHD4-2,IGHD4-3,IGHD4-4,IGHD5-1,\
            IGHD5-2,IGHD5-3,IGHD6-1,IGHD6-2,\
            IGHD6-3,IGHD6-4,IGHD6-5,IGHD6-6,IGHD7-1,none\n')
        out.write(str(cdr_all.d_segs_nc['IGHD1-1']) + ',' + str(cdr_all.d_segs_nc['IGHD1-2']) + ',' + str(cdr_all.d_segs_nc['IGHD1-3']) + ',' + str(cdr_all.d_segs_nc['IGHD1-4']) + ',' + str(cdr_all.d_segs_nc['IGHD1-5'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD1-6']) + ',' + str(cdr_all.d_segs_nc['IGHD1-7']) + ',' + str(cdr_all.d_segs_nc['IGHD1-8']) + ',' + str(cdr_all.d_segs_nc['IGHD2-1']) \
                              + ',' + str(cdr_all.d_segs_nc['IGHD2-2']) + ',' + str(cdr_all.d_segs_nc['IGHD2-3'])  + ',' + str(cdr_all.d_segs_nc['IGHD2-4'])  + ',' + str(cdr_all.d_segs_nc['IGHD2-5'])  + ',' + str(cdr_all.d_segs_nc['IGHD2-6'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD3-1']) + ',' + str(cdr_all.d_segs_nc['IGHD3-2']) + ',' + str(cdr_all.d_segs_nc['IGHD3-3'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD3-4']) + ',' + str(cdr_all.d_segs_nc['IGHD4-1']) + ',' + str(cdr_all.d_segs_nc['IGHD4-2'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD4-3']) + ',' + str(cdr_all.d_segs_nc['IGHD4-4']) + ',' + str(cdr_all.d_segs_nc['IGHD5-1']) + ',' + str(cdr_all.d_segs_nc['IGHD5-2'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD5-3'])+ ',' + str(cdr_all.d_segs_nc['IGHD6-1']) + ',' + str(cdr_all.d_segs_nc['IGHD6-2']) + ',' + str(cdr_all.d_segs_nc['IGHD6-3'])\
                              + ',' + str(cdr_all.d_segs_nc['IGHD6-4']) + ',' + str(cdr_all.d_segs_nc['IGHD6-5'])+ ',' + str(cdr_all.d_segs_nc['IGHD6-6'])+ ',' + str(cdr_all.d_segs_nc['IGHD7-1'])+ ',' + str(cdr_all.d_segs_nc['']) + ','\
                              + ',' + str(cdr_all.d_segs['IGHD1-1']) + ',' + str(cdr_all.d_segs['IGHD1-2']) + ',' + str(cdr_all.d_segs['IGHD1-3']) + ',' + str(cdr_all.d_segs['IGHD1-4']) + ',' + str(cdr_all.d_segs['IGHD1-5'])\
                              + ',' + str(cdr_all.d_segs['IGHD1-6']) + ',' + str(cdr_all.d_segs['IGHD1-7']) + ',' + str(cdr_all.d_segs['IGHD1-8']) + ',' + str(cdr_all.d_segs['IGHD2-1']) \
                              + ',' + str(cdr_all.d_segs['IGHD2-2']) + ',' + str(cdr_all.d_segs['IGHD2-3'])  + ',' + str(cdr_all.d_segs['IGHD2-4'])  + ',' + str(cdr_all.d_segs['IGHD2-5'])  + ',' + str(cdr_all.d_segs['IGHD2-6'])\
                              + ',' + str(cdr_all.d_segs['IGHD3-1']) + ',' + str(cdr_all.d_segs['IGHD3-2']) + ',' + str(cdr_all.d_segs['IGHD3-3'])\
                              + ',' + str(cdr_all.d_segs['IGHD3-4']) + ',' + str(cdr_all.d_segs['IGHD4-1']) + ',' + str(cdr_all.d_segs['IGHD4-2'])\
                              + ',' + str(cdr_all.d_segs['IGHD4-3']) + ',' + str(cdr_all.d_segs['IGHD4-4']) + ',' + str(cdr_all.d_segs['IGHD5-1']) + ',' + str(cdr_all.d_segs['IGHD5-2'])\
                              + ',' + str(cdr_all.d_segs['IGHD5-3'])+ ',' + str(cdr_all.d_segs['IGHD6-1']) + ',' + str(cdr_all.d_segs['IGHD6-2']) + ',' + str(cdr_all.d_segs['IGHD6-3'])\
                              + ',' + str(cdr_all.d_segs['IGHD6-4']) + ',' + str(cdr_all.d_segs['IGHD6-5'])+ ',' + str(cdr_all.d_segs['IGHD6-6'])+ ',' + str(cdr_all.d_segs['IGHD7-1'])+ ',' + str(cdr_all.d_segs[''])) 
    dseg_list_nc = [cdr_all.d_segs_nc['IGHD1-1'], cdr_all.d_segs_nc['IGHD1-2'], cdr_all.d_segs_nc['IGHD1-3'], cdr_all.d_segs_nc['IGHD1-4'], cdr_all.d_segs_nc['IGHD1-5']\
                              , cdr_all.d_segs_nc['IGHD1-6'], cdr_all.d_segs_nc['IGHD1-7'], cdr_all.d_segs_nc['IGHD1-8'], cdr_all.d_segs_nc['IGHD2-1']\
                              , cdr_all.d_segs_nc['IGHD2-2'], cdr_all.d_segs_nc['IGHD2-3'], cdr_all.d_segs_nc['IGHD2-4'], cdr_all.d_segs_nc['IGHD2-5'], cdr_all.d_segs_nc['IGHD2-6']\
                              , cdr_all.d_segs_nc['IGHD3-1'], cdr_all.d_segs_nc['IGHD3-2'], cdr_all.d_segs_nc['IGHD3-3']\
                              , cdr_all.d_segs_nc['IGHD3-4'], cdr_all.d_segs_nc['IGHD4-1'], cdr_all.d_segs_nc['IGHD4-2']\
                              , cdr_all.d_segs_nc['IGHD4-3'], cdr_all.d_segs_nc['IGHD4-4'], cdr_all.d_segs_nc['IGHD5-1'], cdr_all.d_segs_nc['IGHD5-2']\
                              , cdr_all.d_segs_nc['IGHD5-3'], cdr_all.d_segs_nc['IGHD6-1'], cdr_all.d_segs_nc['IGHD6-2'], cdr_all.d_segs_nc['IGHD6-3']\
                              , cdr_all.d_segs_nc['IGHD6-4'], cdr_all.d_segs_nc['IGHD6-5'], cdr_all.d_segs_nc['IGHD6-6'], cdr_all.d_segs_nc['IGHD7-1'], cdr_all.d_segs_nc['']]
    dseg_list = [cdr_all.d_segs['IGHD1-1'], cdr_all.d_segs['IGHD1-2'], cdr_all.d_segs['IGHD1-3'], cdr_all.d_segs['IGHD1-4'], cdr_all.d_segs['IGHD1-5']\
                              , cdr_all.d_segs['IGHD1-6'], cdr_all.d_segs['IGHD1-7'], cdr_all.d_segs['IGHD1-8'], cdr_all.d_segs['IGHD2-1'] \
                              , cdr_all.d_segs['IGHD2-2'], cdr_all.d_segs['IGHD2-3'], cdr_all.d_segs['IGHD2-4'], cdr_all.d_segs['IGHD2-5'], cdr_all.d_segs['IGHD2-6']\
                              , cdr_all.d_segs['IGHD3-1'], cdr_all.d_segs['IGHD3-2'], cdr_all.d_segs['IGHD3-3']\
                              , cdr_all.d_segs['IGHD3-4'], cdr_all.d_segs['IGHD4-1'], cdr_all.d_segs['IGHD4-2']\
                              , cdr_all.d_segs['IGHD4-3'], cdr_all.d_segs['IGHD4-4'], cdr_all.d_segs['IGHD5-1'], cdr_all.d_segs['IGHD5-2']\
                              , cdr_all.d_segs['IGHD5-3'], cdr_all.d_segs['IGHD6-1'], cdr_all.d_segs['IGHD6-2'], cdr_all.d_segs['IGHD6-3']\
                              , cdr_all.d_segs['IGHD6-4'], cdr_all.d_segs['IGHD6-5'], cdr_all.d_segs['IGHD6-6'], cdr_all.d_segs['IGHD7-1'], cdr_all.d_segs['']]
    return [dseg_list_nc, dseg_list]



def plot_dsegs(infile, outfile, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    data = analyse_dsegs(folders, outfile, regex, frequency)
    ind = np.arange(len(data[0]))
    labels = ['D1-1','D1-2','D1-3','D1-4','D1-5','D1-6','D1-7','D1-8',\
            'D2-1','D2-2','D2-3','D2-4','D2-5','D2-6',\
            'D3-1','D3-2','D3-3','D3-4','D4-1',\
            'D4-2','D4-3','D4-4','D5-1',\
            'D5-2','D5-3','D6-1','D6-2',\
            'D6-3','D6-4','D6-5','D6-6','D7-1','none']
    font = fm.FontProperties(family = 'Sans-Serif', size = 20)
    font2 = fm.FontProperties(family = 'Sans-Serif', size = 25)
    font1 = fm.FontProperties(family = 'Sans-Serif', size = 20, weight = 'bold')
    values = range(9)
    jet = cm = plt.get_cmap('CMRmap') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    fig = plt.figure(facecolor = 'white')
    fig.set_size_inches(10, 10)
    colorVal = scalarMap.to_rgba(values[4])
    colorVal2 = scalarMap.to_rgba(values[2])
    plt.bar(ind, data[0], color = colorVal, label = 'non-canonical')
    plt.bar(ind, data[1], color = colorVal2, bottom = data[0], label = 'canonincal')
    plt.ylim([0,0.8])
    plt.xlim([0,33])
    plt.xticks(ind, labels, rotation=90, rotation_mode="anchor")
    plt.tick_params(axis = 'x', pad = 40, which = 'both', bottom = 'off', top = 'off')
    ax = plt.gca()
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontproperties(font)
    xticks = ax.xaxis.get_major_ticks()
    for t in xticks[8:14]:
        t.label1.set_fontproperties(font1)
    plt.tight_layout()
    plt.savefig(outfile + '.png')



                                          
                
                        





