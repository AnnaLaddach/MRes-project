#Script to analyse and plot the relative frequency of cysteine spacings in CDR3s containing two or more cysteines.
#A textfile listing IMGT output folders and output file name are taken as input. 
#Optionally sequence frequency information can be included - a regex must be supplied to extract this from the sequence ID.


import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import helper_functions
import matplotlib.font_manager as fm
import matplotlib.colors as colors
import matplotlib.cm as cmx


def main(argv):
    if len(argv) == 3:
        plot_space(argv[1], argv[2])
    elif len(argv) == 5:
	plot_space(argv[1], argv[2], argv[3], argv[4])
    else:
        print 'usage: cysteine_spacing.py IMGT_folder_list outfile '
	print 'optional arguments: frequency_regex frequency(TRUE/FALSE)'
        sys.exit()


def space(folders, outfile, regex = None, frequency = False):
    space02_all = {}
    space3_all = {}
    space4_all = {}
    space56_all = {}
    space_large_all = {}
    for i in range(1,51):
        space02_all[i] = []
        space3_all[i] = []
        space4_all[i] = []
        space56_all[i] = []
        space_large_all[i] = []
    for folder in folders:
        with open(folder[0] + '/5_AA-sequences.txt') as f:
            space02 = {}
            space3 = {}
            space4 = {}
            space56 = {}
            space_large ={}
            lengths = {}
            for i in range(1,51):
                space02[i] = 0
                space3[i] = 0
                space4[i] = 0
                space56[i] = 0
                space_large[i] = 0
                lengths[i] = 0
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                if row['Functionality'].find('productive') == 0:
                    try:
                        if frequency:
                            pat = re.compile(regex)
                            info = pat.match(row['Sequence ID'])
                            freq = int(info.group(1))
                        else:
                            freq = 1
                        if row['CDR3-IMGT'].count('C') > 1:
                            length = len(row['CDR3-IMGT'])
                            space = row['CDR3-IMGT'].rfind('C') - row['CDR3-IMGT'].find('C') - 1
                            if space < 3:
                                space02[length] += freq
                            elif space == 3:
                                space3[length] += freq
                            elif space == 4:
                                space4[length] += freq
                            elif space in [5,6]:
                                space56[length] += freq
                            else:
                                space_large[length] += freq
                            lengths[length] += freq
                    except:
                        pass
            for i in range(1,51):
                try:
                    space02_all[i].append(space02[i]/float(lengths[i]))
                    space3_all[i].append(space3[i]/float(lengths[i]))
                    space4_all[i].append(space4[i]/float(lengths[i]))
                    space56_all[i].append(space56[i]/float(lengths[i]))
                    space_large_all[i].append(space_large[i]/float(lengths[i]))
                except:
                    pass
    data = []
    with open(outfile + '.txt', 'w') as out:
        out.write('CDR3_length,space_0-2,space_3,space4,space_5-6,space_7+\n')
        for i in range(1,51):
            out.write(str(i) + ',' +  str(np.mean(space02_all[i]))  + ',' +  str(np.mean(space3_all[i]))  + ','\
                      +  str(np.mean(space4_all[i]))  + ',' +  str(np.mean(space56_all[i]))  + ',' +  str(np.mean(space_large_all[i])) + '\n')
            data.append([i,np.mean(space02_all[i]),np.mean(space3_all[i]),np.mean(space4_all[i]),np.mean(space56_all[i]),np.mean(space_large_all[i])])
    return data



            
def plot_space(infile, outfile, regex = None, frequency = False):
    if frequency:
        outfile = outfile + '_freq'
    folders = helper_functions.create_list(infile)
    data = space(folders, outfile, regex, frequency)
    arr = np.array(data)
    fig = plt.figure(facecolor = 'white')
    values = range(6)
    c = cm = plt.get_cmap('CMRmap') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=c)
    plt.ylabel('relative frequency of cysteine spacing')
    plt.xlabel('CDR3 length (AA)')
    plt.plot(arr[:,0], arr[:,1], color =  scalarMap.to_rgba(values[0]), label = '0-2',linewidth=1.5)
    plt.scatter(arr[:,0], arr[:,1], color =  scalarMap.to_rgba(values[0]))
    plt.plot(arr[:,0], arr[:,2], color =  scalarMap.to_rgba(values[1]), label = '3', linewidth=1.5)
    plt.scatter(arr[:,0], arr[:,2], color =  scalarMap.to_rgba(values[1]))
    plt.plot(arr[:,0], arr[:,3], color =  scalarMap.to_rgba(values[2]), label = '4', linewidth=1.5)
    plt.scatter(arr[:,0], arr[:,3], color =  scalarMap.to_rgba(values[2]))
    plt.plot(arr[:,0], arr[:,4], color =  scalarMap.to_rgba(values[3]), label = '5-6', linewidth=1.5)
    plt.scatter(arr[:,0], arr[:,4], color =  scalarMap.to_rgba(values[3]))
    plt.plot(arr[:,0], arr[:,5], color =  scalarMap.to_rgba(values[4]), label = '>7', linewidth=1.5)
    plt.scatter(arr[:,0], arr[:,5], color =  scalarMap.to_rgba(values[4]))
    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlim(10,31)
    plt.ylim([0,1])
    plt.savefig(outfile + '.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

if __name__=="__main__":
     main(sys.argv) 



                          
                    

                            
                            
                            
                    

                            
                            
                            
                    

                            
                            
                            
                    
