#Script to analyse and plot amino acid substitutions from germline IGHD2-2*, IGHD2-15*, IGHD2-8*01, IGHD2-8*02, IGHD3-3*, IGHD3-9*, IGHD3-10* and IGHD3-22* gene segments
#A textfile listing IMGT output folders, output file name and germline gene segment name are taken as input. 
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID.

import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import helper_functions
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.font_manager as fm

def main(argv):
    if len(argv) == 4:
        plot_substitutions(argv[1], argv[2], argv[3])
    elif len(argv) == 5:
	plot_substitutions(argv[1], argv[2], argv[3], argv[4], argv[5])
    else:
        print 'usage: analyse_substitutions.py IMGT_folder_list outfile germline_d_segment'
	print 'optional arguments: frequency_regex frequency(TRUE/FALSE)'
        sys.exit()


#global constants
seq_dict = {'IGHD2-2*':'SSTS','IGHD2-15*':'SGGS', 'IGHD2-8*01':'TNGV', 'IGHD2-8*02':'TNGV',\
            'IGHD3-3*':'FWSG', 'IGHD3-9*':'ILTG', 'IGHD3-10*':'GSGS', 'IGHD3-22*':'DSSG'}
amino_acids = ['P','G','A','V', 'L', 'I', 'F', 'Y', 'W',\
        'M', 'C', 'S', 'T', 'N', 'Q', 'H', 'K', 'R', 'D', 'E']

class Space:
    def __init__(self):
        self.pos = [0,0,0,0]
        self.mutated = [0,0,0,0]
        self.aminoacids = []
        for i in range(4):
            self.aminoacids.append({})

    def increment(self, index, frequency):
        self.pos[index] += frequency

    def add_amino_acid(self, amino_acid, position, frequency):
        if amino_acid in self.aminoacids[position].keys():
            self.aminoacids[position][amino_acid] += frequency
        else:
            self.aminoacids[position][amino_acid] = frequency
        self.mutated[position] += frequency

class Space_all:
    def __init__(self):
        self.mutated = []
        for i in range(4):
            self.mutated.append([])
        self.aminoacids = []
        for i in range(4):
            self.aminoacids.append({})
        self.aminoacids_means = []
        for i in range(4):
            self.aminoacids_means.append({})
        self.aminoacids_sd = []
        for i in range(4):
            self.aminoacids_sd.append({})
        self.mutated_means = []
        self.mutated_sd = []
        
    def increment(self, index, mutated):
        self.mutated[index].append(mutated)

    def add_amino_acid(self, amino_acid, position, frequency):
        if amino_acid in self.aminoacids[position].keys():
            self.aminoacids[position][amino_acid].append(frequency)
        else:
            self.aminoacids[position][amino_acid] = [frequency]

    def calculate_stats(self):
        i = 0
        for aa_dict in self.aminoacids:
            for aa in aa_dict.keys():
                while len(aa_dict[aa]) < len(self.mutated[0]):
                    aa_dict[aa].append(0)
                self.aminoacids_means[i][aa] = np.mean(aa_dict[aa])
                self.aminoacids_sd[i][aa] = np.std(aa_dict[aa])
            i += 1
        for ind in self.mutated:
            self.mutated_means.append(np.mean(ind))
            self.mutated_sd.append(np.std(ind)) 
                

def substitutions(folders, outfile, germline, regex = None, frequency = False):
    space_all = Space_all()
    n = 0
    try:
        seq = seq_dict[germline]
    except:
        print 'Make sure germline is one of IGHD2-2*, IGHD2-15*, IGHD2-8*01, IGHD2-8*02, IGHD3-3*, IGHD3-9*, IGHD3-10*, IGHD3-22*'
        return False
    if germline in ['IGHD2-2*','IGHD2-15*','IGHD2-8*01','IGHD2-8*02']:
        pat = re.compile('(?:tgt|tgc)(.{12})(?:tgt|tgc)')
    elif germline in ['IGHD3-3*','IGHD3-9*']:
        pat = re.compile('(?:gat|gac)(.{12})(?:tat|tac)')
    elif germline in ['IGHD3-10*', 'IGHD3-22*']:
        pat = re.compile('(?:tat|tac)(.{12})(?:tat|tac)')
    else:
        print 'Make sure you have germline is one of IGHD2-2*, IGHD2-15*, IGHD2-8*01, IGHD2-8*02, IGHD3-3*, IGHD3-9*, IGHD3-10*, IGHD3-22*'
        return False
    for folder in folders:
        with open(folder[0] + '/6_Junction.txt') as f:
            space = Space()
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'].find('productive') == 0:
                    if row['D-GENE and allele'].find(germline)!= -1 and row['D-REGION reading frame'] == '2':
                        try:
                            info = pat.search(row['D-REGION'])
                            nucs = info.group(1)
                            seq1 = Seq(nucs, generic_dna)
                            letters = str(seq1.translate())
                            if frequency:
                                pat = re.compile(regex)
                                info = pat.match(row['Sequence ID'])
                                freq = int(info.group(1))
                            else:
                                freq = 1
                            n += freq
                            i = 0
                            for letter in letters:
                                if letter != seq[i]:
                                    space.add_amino_acid(letters[i], i, freq)
                                space.increment(i, freq)
                                i += 1
                        except:
                            pass
                  
            i = 0
            for aa_dict in space.aminoacids:
                for aa in aa_dict.keys():
                    space_all.add_amino_acid(aa, i, aa_dict[aa]/float(space.mutated[i]))
                try:
                    space_all.increment(i, space.mutated[i]/float(space.pos[i]))
                except:
                    space_all.increment(i, 0)
                i += 1
    space_all.calculate_stats()
    data = [[],[],[],[]]
    with open(outfile + '_mutations_' + str(n) + '.csv', 'w') as out:
        with open(outfile + '_' + str(n) + '.csv', 'w') as out1:
            out.write('position,total_fraction_mutated,standard_dev_fraction_mutated,' + ','.join(amino_acids) + ',sdev ' + ',sdev '.join(amino_acids) + '\n')
            i = 1
            for aa_dict in space_all.aminoacids:
                out.write(str(i) +',')
                out.write(str(space_all.mutated_means[i-1]) + ',')
                data[i-1].append(space_all.mutated_means[i-1])
                out.write(str(space_all.mutated_sd[i-1]))
                data[i-1].append(space_all.mutated_sd[i-1])
                out1.write(','.join(str(x) for x in space_all.mutated[i-1]) +'\n')
                for aa in amino_acids:
                    if aa in aa_dict.keys():
                        out.write(',' + str(space_all.aminoacids_means[i-1][aa]))
                        data[i-1].append(space_all.aminoacids_means[i-1][aa])
                    else:
                        out.write(',0')
                        data[i-1].append(0)
                for aa in amino_acids:
                    if aa in aa_dict.keys():
                        out.write(',' + str(space_all.aminoacids_sd[i-1][aa]))
                        data[i-1].append(space_all.aminoacids_sd[i-1][aa])
                    else:
                        out.write(',0')
                        data[i-1].append(0)
                i += 1
                out.write('\n')
    return (data, n)


def plot_substitutions(infile, outfile, germline, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    amino_acids = ['P','G','A','V', 'L', 'I', 'F', 'Y', 'W',\
        'M', 'C', 'S', 'T', 'N', 'Q', 'H', 'K', 'R', 'D', 'E']
    seq = seq_dict[germline]
    folders = helper_functions.create_list(infile)
    data = substitutions(folders, outfile, germline, regex, frequency)
    if not data:
        return
    arr = np.array(data[0])
    ind = np.arange(4)
    values = range(4)
    c = cm = plt.get_cmap('CMRmap') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=c)
    font = fm.FontProperties(family = 'Sans-Serif', size = 25)
    font1 = fm.FontProperties(family = 'Sans-Serif', size = 15)
    font2 = fm.FontProperties(family = 'Sans-Serif', size = 50)

    prob_list_d = [1.509482E-01, 2.112400E+00, 1.241095E+00, 4.076074E-01, 3.020502E-02, 3.916999E-02, 1.031734E-01,\
                   9.078338E-01, 1.423897E-02, 3.357577E-02, 2.340019E-02, 5.822988E-01, 1.516881E-01, 7.185182E+00,\
                   1.843856E-01, 1.126995E+01, 6.622772E-02, 5.795374E-02, 0, 7.426619E+00]
    total_d = sum(prob_list_d)
    d_ab = [prob/float(total_d) for prob in prob_list_d]
    prob_list_f = [3.556198E-01, 1.984772E-01, 1.130146E-01, 1.643946E+00, 3.390618E+00, 2.714705E+00, 0, 7.074464E+00,\
                   2.199856E-01, 3.947940E+00, 1.993818E-01, 9.813064E-01, 3.674486E-02, 1.955446E-01, 1.496610E-03,\
                   5.642309E+00, 4.649035E-03, 8.208677E-02, 1.031734E-01, 5.288625E-02]
    total_f = sum(prob_list_f)
    f_ab = [prob/float(total_f) for prob in prob_list_f]
    prob_list_g = [7.915056E-02, 0, 1.851951E+00, 4.228200E-01, 2.157936E-02, 3.560482E-02, 1.984772E-01,\
                   2.295842E-01, 1.386853E-01, 2.187264E-02, 5.891123E-02, 9.984215E-01, 7.862767E-02,\
                   4.868901E-01, 5.516340E-02, 9.139518E+00, 5.504400E-02, 1.089400E+00, 2.112400E+00, 1.389370E+00]
    total_g = sum(prob_list_g)
    g_ab = [prob/float(total_g) for prob in prob_list_g]
    prob_list_i = [5.016568E-02, 3.560482E-02, 1.140412E-01, 7.926675E+00, 1.601123E+00, 0, 2.714705E+00, 2.903165E-01, 6.318748E-02,\
                   1.096748E+01, 6.594967E-04, 6.007607E-01, 3.087152E+00, 1.762721E+00, 6.712736E-06, 4.706586E+00, 5.059793E-01,\
                   3.245175E-01, 3.916999E-02, 1.029959E-04]
    total_i = sum(prob_list_i)
    i_ab = [prob/float(total_i) for prob in prob_list_i]
    prob_list_l = [1.149931E+00, 2.157936E-02, 6.969101E-02, 3.595310E+00, 0, 1.601123E+00, 3.390618E+00, 1.521320E-01, 3.378544E-01,\
                   5.647985E+00, 6.079219E-03, 1.580539E-01, 5.702792E-01, 2.769442E-02, 6.802781E-01, 5.879103E+00, 2.158451E-02, \
                   3.932002E-01, 3.020502E-02, 1.283121E-03]
    total_l = sum(prob_list_l)
    l_ab = [prob/float(total_l) for prob in prob_list_l]
    prob_list_n = [7.342554E-03, 4.868901E-01, 9.291290E-02, 4.931206E-02, 2.769442E-02, 1.762721E+00, 1.955446E-01,\
                   2.622763E+00, 7.310937E-04, 3.972173E-02, 1.374268E-06, 7.348346E+00, 2.147175E+00, 0, 8.705989E-02,\
                   1.405444E+01, 6.190318E+00, 7.829130E-01, 7.185182E+00, 5.038373E-01]
    total_n = sum(prob_list_n)
    n_ab = [prob/float(total_n) for prob in prob_list_n]
    prob_list_s = [1.284651E+00, 9.984215E-01, 9.988358E-01, 7.477041E-02, 1.580539E-01, 6.007607E-01, 9.813064E-01,\
                8.104751E-01, 1.385142E-01, 1.861354E-02, 2.639482E-01, 0, 3.058575E+00, 7.348346E+00, 5.906405E-04,\
                5.439116E+00, 8.688405E-02, 1.926435E+00, 5.822988E-01, 6.776709E-02]
    total_s = sum(prob_list_s)
    s_ab = [prob/float(total_s) for prob in prob_list_s]
    prob_list_t = [9.057112E-01, 7.862767E-02, 2.912317E+00, 2.166054E-01, 5.702792E-01, 3.087152E+00, 3.674486E-02,\
                   9.984255E-02, 1.412361E-02, 1.415612E+00, 3.225214E-06, 3.058575E+00, 0, 2.147175E+00, 1.202094E-01,\
                   3.443285E+00, 1.039298E+00, 1.135258E+00, 1.516881E-01, 6.016624E-02]
    total_t = sum(prob_list_t)
    t_ab = [prob/float(total_t) for prob in prob_list_t]
    prob_list_v = [2.217442E-01, 4.228200E-01, 3.774477E+00, 0, 3.595310E+00, 7.926675E+00, 1.643946E+00, 5.010635E-01,\
                   9.663569E-02, 4.396720E+00, 2.243512E-02, 7.477041E-02, 2.166054E-01, 4.931206E-02, 9.047737E-03,\
                   6.890244E+00, 3.493440E-02, 1.366145E-01, 4.076074E-01, 5.795409E-01]
    total_v = sum(prob_list_v)
    v_ab = [prob/float(total_v) for prob in prob_list_v]
    prob_list_w = [5.516074E-03, 1.386853E-01, 7.939549E-02, 9.663569E-02, 3.378544E-01, 6.318748E-02, 2.199856E-01,\
                   6.121284E-01, 0, 1.011149E-01, 4.440833E-01, 1.385142E-01, 1.412361E-02, 7.310937E-04, 4.332983E-05,\
                   7.013890E+00, 8.024263E-03, 5.724286E-01, 1.423897E-02, 2.252612E-02]
    total_w = sum(prob_list_w)
    w_ab = [prob/float(total_w) for prob in prob_list_w]

    prob_list_d_lg = [0.394456, 0.844926, 0.395144, 0.037967, 0.015076, 0.010690, 0.017416, 0.135107, 0.029890,\
                      0.025548, 0.062556, 1.240275, 0.425860, 5.076149, 0.523386, 0.927114, 0.282959, 0.123954, 0, 5.243870]    
    total_d_lg = sum(prob_list_d_lg)
    d_lg = [prob/float(total_d_lg) for prob in prob_list_d_lg]
    prob_list_f_lg = [0.094464, 0.089586, 0.253701, 0.654683, 2.592692, 1.112727, 0, 7.803902, 2.457121, 1.798853,\
                      1.105251, 0.361819, 0.165001, 0.089525, 0.035855, 0.682139, 0.023918, 0.052722, 0.017416, 0.018811]
    total_f_lg = sum(prob_list_f_lg)
    f_lg = [prob/float(total_f_lg) for prob in prob_list_f_lg]
    prob_list_g_lg = [0.196961, 0, 2.066040, 0.076701, 0.044261, 0.008705, 0.089586, 0.054679, 0.268491, 0.139538,\
                      0.569265, 1.739990, 0.129836, 1.437645, 0.267959, 0.311484, 0.296636, 0.390192, 0.844926, 0.348847]
    total_g_lg = sum(prob_list_g)
    g_lg = [prob/float(total_g_lg) for prob in prob_list_g_lg]
    prob_list_i_lg = [0.078281, 0.008705, 0.149830, 10.649107, 4.145067, 0, 1.112727, 0.232523, 0.111660, 4.273607, 0.320627,\
                      0.064105, 1.033739, 0.191503, 0.072854, 0.108882, 0.159069, 0.126991, 0.010690, 0.044265]
    total_i_lg = sum(prob_list_i_lg)
    i_lg = [prob/float(total_i_lg) for prob in prob_list_i_lg]
    prob_list_l_lg = [0.249060, 0.044261, 0.395337, 1.702745, 0, 4.145067, 2.592692, 0.299648, 0.619632, 6.312358, 0.594007,\
                      0.182287, 0.302936, 0.068427, 0.582457, 0.366317, 0.137500, 0.301848, 0.015076, 0.069673]
    total_l_lg = sum(prob_list_l_lg)
    l_lg = [prob/float(total_l_lg) for prob in prob_list_l_lg]
    prob_list_n_lg = [0.161787, 1.437645, 0.276818, 0.083688, 0.068427, 0.191503, 0.089525, 0.612025, 0.045376,\
                      0.371004, 0.528768, 4.008358, 2.000679, 0, 1.695752, 4.509238, 2.145078, 0.751878, 5.076149, 0.541712]
    total_n_lg = sum(prob_list_n_lg)
    n_lg = [prob/float(total_n_lg) for prob in prob_list_n_lg]
    prob_list_s_lg = [1.338132, 1.739990, 4.727182, 0.098369, 0.182287, 0.064105, 0.361819, 0.400547, 0.248862,\
                      0.346960, 2.784478, 0, 6.472279, 4.008358, 1.223828, 0.990012, 0.748683, 0.858151, 1.240275, 0.611973]
    total_s_lg = sum(prob_list_s_lg)
    s_lg = [prob/float(total_s_lg) for prob in prob_list_s_lg]
    prob_list_t_lg = [0.571468, 0.129836, 2.139501, 2.188158, 0.302936, 1.033739, 0.165001, 0.245841, 0.140825,\
                      2.020366, 1.143480, 6.472279, 0, 2.000679, 1.080136, 0.584262, 1.136863, 0.578987, 0.425860, 0.604545]
    total_t_lg = sum(prob_list_t_lg)
    t_lg = [prob/float(total_t_lg) for prob in prob_list_t_lg]
    prob_list_v_lg = [0.296501, 0.076701, 2.547870, 0, 1.702745, 10.649107, 0.654683, 0.249313, 0.189510,\
                      1.898718, 1.959291, 0.098369, 2.188158, 0.083688, 0.210332, 0.119013, 0.185202, 0.170887,\
                      0.037967, 0.245034]
    total_v_lg = sum(prob_list_v_lg)
    v_lg = [prob/float(total_v_lg) for prob in prob_list_v_lg]
    prob_list_w_lg = [0.095131, 0.268491, 0.180717, 0.189510, 0.619632, 0.111660, 2.457121, 3.151815, 0, 0.696175,\
                      0.670128, 0.248862, 0.140825, 0.045376, 0.236199, 0.597054, 0.049906, 0.593607, 0.029890, 0.077852]
    total_w_lg = sum(prob_list_w_lg)
    w_lg = [prob/float(total_w_lg) for prob in prob_list_w_lg]

    
    ab = {'D': d_ab, 'F' : f_ab, 'G' : g_ab, 'I': i_ab, 'L': l_ab, 'N': n_ab, 'S': s_ab, 'T': t_ab, 'V': v_ab, 'W': w_ab}
    lg = {'D': d_lg, 'F' : f_lg, 'G' : g_lg, 'I': i_lg, 'L': l_lg, 'N': n_lg, 'S': s_lg, 'T': t_lg, 'V': v_lg, 'W': w_lg}
    fig = plt.figure(facecolor = 'white')
    plt.bar(ind, arr[:,0], yerr = arr[:,1], align = 'center')
    labels = [1,2,3,4]
    width = 0.3
    plt.xlabel('position')
    plt.ylabel('relative frequency of mutations')
    plt.xticks(ind, labels)
    plt.ylim([0,0.6])
    plt.savefig(outfile + '_' + str(data[1]) + '_mutations.png')
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, sharex=True, sharey=False, facecolor = 'white', figsize = (10,10))
    ind2 = np.arange(20)
    plt.setp((ax0, ax1, ax2, ax3), xticks = ind2, xticklabels = amino_acids) 
    c1 = color = scalarMap.to_rgba(values[0])
    c2 = color = scalarMap.to_rgba(values[1])
    c3 = color = scalarMap.to_rgba(values[2])
    ax0.bar(ind2 + 0.05, arr[0,2:22], width, color = c3, label = 'observed substitutions')
    ax0.bar(ind2 + 0.35, ab[seq[0]], width, color = c2, label = 'AB model')
    ax0.bar(ind2 + 0.65, lg[seq[0]], width,color = c1, label = 'LG model')
    ax0.set_ylim([0,0.6])
    ax1.bar(ind2 + 0.05, arr[1,2:22], width, color = c3)
    ax1.bar(ind2 + 0.35, ab[seq[1]], width, color = c2)
    ax1.bar(ind2 + 0.65, lg[seq[1]], width, color = c1)
    ax1.set_ylim([0,0.6])
    ax2.bar(ind2 + 0.05, arr[2,2:22], width, color = c3)
    ax2.bar(ind2 + 0.35, ab[seq[2]], width, color = c2)
    ax2.bar(ind2 + 0.65, lg[seq[2]], width, color = c1)
    ax2.set_ylim([0,0.6])
    ax3.bar(ind2 +0.05, arr[3,2:22], width, color = c3)
    ax3.bar(ind2 + 0.35, ab[seq[3]], width, color = c2)
    ax3.bar(ind2 + 0.65, lg[seq[3]], width, color = c1)
    ax3.set_ylim([0,0.6])
    for label in (ax3.get_xticklabels()):
        label.set_fontproperties(font)
    for label in (ax0.get_yticklabels()):
        label.set_fontproperties(font1)
    for label in (ax1.get_yticklabels()):
        label.set_fontproperties(font1)
    for label in (ax2.get_yticklabels()):
        label.set_fontproperties(font1)
    for label in (ax3.get_yticklabels()):
        label.set_fontproperties(font1)
    ax0.set_ylabel('1', rotation = 0, fontproperties = font2)
    ax1.set_ylabel('2', rotation = 0, fontproperties = font2)
    ax2.set_ylabel('3', rotation = 0, fontproperties = font2)
    ax3.set_ylabel('4', rotation = 0, fontproperties = font2)
    ax0.yaxis.labelpad = 30
    ax1.yaxis.labelpad = 30
    ax2.yaxis.labelpad = 30
    ax3.yaxis.labelpad = 30
    ax0.set_xticks(ind2, minor=True)
    ax0.grid(True, which = 'minor')
    ax1.set_xticks(ind2, minor=True)
    ax1.grid(True, which = 'minor') 
    ax2.set_xticks(ind2, minor=True)
    ax2.grid(True, which = 'minor')
    ax3.set_xticks(ind2, minor=True)
    ax3.grid(True, which = 'minor')
    lgd = ax0.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop = font)
    plt.xticks(ind2 + 0.5, amino_acids)
    plt.tight_layout()
    plt.savefig(outfile + '_' + str(data[1]) + '_amino_acids.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

    
if __name__=="__main__":
     main(sys.argv) 
 
                            

                        

                                                
                    





            

    
