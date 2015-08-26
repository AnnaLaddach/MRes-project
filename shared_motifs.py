#Script to enumerate and plot the number of shared motifs between individuals arising between two conserved cysteines of IGHD2 germline gene segments.
#A textfile listing IMGT output folders and output file, germline gene segment and frequency threshold for accepting a motif are taken as input. 
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID.
#N.B. using frequency information will only affect the output if a threshold greater than 0 is used.

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
    if len(argv) == 5:
        run_shared_motifs(argv[1], argv[2], argv[3], int(argv[4]))
    elif len(argv) == 7:
	run_shared_motifs(argv[1], argv[2], argv[3], int(argv[4]), argv[5], argv[6])
    else:
        print 'usage: enumerate_motifs.py IMGT_folder_list outfile germline_gene_segment threshold'
	print 'germline gene segment can be one of \'IGHD2-8*01\', \'IGHD2-8*02\', \'IGHD2-2*\', \'IGHD2-15*\' or \'IGHD2-21*'
	print 'optional arguments: frequency_regex frequency(TRUE/FALSE)'
        sys.exit()

def shared_motifs(folders, outfile, germline, threshold, regex = None, frequency = False):
    if germline in ['IGHD2-8*01', 'IGHD2-8*02', 'IGHD2-2*', 'IGHD2-15*']:
        pat = re.compile('C(....)C')
    elif germline == 'IGHD2-21*':
        pat = re.compile('C(...)C')
    all_seqs = {}
    for folder in folders:
        with open(folder[0] + '/5_AA-sequences.txt') as f:
            patterns = set()
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'].find('productive') == 0:
                    if row['D-GENE and allele'].find(germline) != -1:
                        if frequency:
                            p = re.compile(regex)
                            info1 = p.match(row['Sequence ID'])
                            freq = int(info1.group(1))
                        else:
                            freq = 1
                        if freq > threshold:
                            try:
                                info = pat.search(row['CDR3-IMGT'])
                                patterns.add(info.group(1))
                            except:
                                pass
            for motif in patterns:
                if motif in all_seqs.keys():
                    all_seqs[motif] += freq
                else:
                    all_seqs[motif] = freq
    proportions = {}
    for i in range(1, len(folders) + 1):
        proportions[i] = []
        for key in all_seqs.keys():
            if all_seqs[key] == i:
                proportions[i].append(key)
    data = []
    with open(outfile + '.csv', 'w') as out:
        out.write('number of individuals,fraction,number of motifs,cumulative_number,motifs\n')
        cumulative = 0
        for i in range(len(folders), 0, -1):
            cumulative += len(proportions[i])
            out.write(str(i) + ',' + str(i/float(len(folders))) + ',' \
                      + str(len(proportions[i])) + ',' + str(cumulative) + ',' + ','.join(proportions[i]) + '\n')
            data.append([i,cumulative])
    return data
            


def run_shared_motifs(infile, outfile, germline, threshold, regex = None, frequency = False):
        if frequency:
            outfile = outfile + '_freq'
        folders = helper_functions.create_list(infile)
        font = fm.FontProperties(family = 'Sans-Serif', size = 20)
        font2 = fm.FontProperties(family = 'Sans-Serif', size = 25)
        data = shared_motifs(folders, outfile, germline, threshold, regex, frequency)
        arr = np.array(data)
        values = range(9)

        c = cm = plt.get_cmap('CMRmap') 
        cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=c)
        fig = plt.figure(facecolor = 'white')
        fig.set_size_inches(10, 8)
        colorVal = scalarMap.to_rgba(values[0]) 
        plt.plot(arr[:,0],arr[:,1], color = colorVal, linewidth = 2.0)
        plt.scatter(arr[:,0],arr[:,1], color = colorVal)
        ax = plt.gca()
        ax.set_xlim(0,29)
        ax.set_ylim(0,2600)
        ax.set_xlabel('number of individuals', fontproperties = font2)
        ax.set_ylabel('number of common motifs',fontproperties = font2)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontproperties(font)
        yticks = ax.yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)
        
        plt.savefig(outfile + '.png')

if __name__=="__main__":
     main(sys.argv) 
            

