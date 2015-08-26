#pipeline for the creation of lineage trees

import re
import os
import sys
import glob
import subprocess
import tree_functions3
from ete2 import Tree, TreeStyle, TextFace,  NodeStyle, faces, AttrFace, CircleFace

def main(argv):
    if len(argv) != 3:
        print 'a list of cluster names and cluster_seqs file must be specified'

    pat = re.compile('(/.*/)([^/]+)CDR\._cluster_.*?.txt')
    file_process = pat.match(argv[2])
    folder = file_process.group(1)
    file_name = file_process.group(2)
    clusters = argv[1].split(';')
 
    with open(argv[2]) as cluster_seqs:
        cluster_seqs = cluster_seqs.read()

    with open(folder + file_name + '.fasta') as seqs_file:
        seqs = seqs_file.read()
            
    for cluster in clusters:
	#find V reqions of clustered CDR3s and extract first timepoint
        if not glob.glob(folder + cluster + '/all_seqs.fasta'):
            os.system('mkdir ' + folder + cluster)
	    os.system('mkdir ' + folder + cluster + '/changes')
            tree_functions3.write_fasta_files(cluster, cluster_seqs, seqs, folder + cluster + '/first_timepoint.fasta',\
                                                folder + cluster + '/all_seqs.fasta')

	#infer UCA
        if glob.glob(folder + cluster + '/nucs.csv') and not glob.glob(folder + cluster + '/germlines.fasta'):
            os.system('python /home/anna/Documents/software/for_anna/GermlineFromIMGT.py ' + folder + cluster + '/nucs.csv /home/anna/Documents/MRes/germline/imgt_germlines.fasta "Homo sapiens" ' +  folder + cluster + '/germlines.fasta cj')
    
    for cluster in clusters:        
        if glob.glob(folder + cluster + '/germlines.fasta') and not glob.glob(folder + cluster + '/all_seqs.nt_ali.fasta'):
            ready = ''
            while ready != 'Y':
                ready = raw_input('Please enter Y if the germline sequence has been added to all_seqs.fasta, else enter N: ')
                if ready == 'N':
                    sys.exit()

            if ready =='Y':
		#perform a codon based alignment of family members
                os.system('perl /home/anna/Documents/MRes/software/translatorx_vLocal.pl -i ' + folder + cluster + '/all_seqs.fasta -o ' + folder + cluster + '/all_seqs -p M -c 1')

    for cluster in clusters:        
        if glob.glob(folder + cluster + '/all_seqs.nt_ali.fasta') and not glob.glob(folder + cluster + '/tree.fasta.treefile'):
            ready1 = ''
            while ready1 != 'Y':
                ready1 = raw_input('Please enter Y if the alignment files have been manually curated, else enter N: ')
                if ready1 == 'N':
                    sys.exit()

            if ready1 == 'Y':
                clone_info = tree_functions3.process_fasta_for_tree(folder + cluster + '/all_seqs.nt_ali.fasta', folder + cluster + '/tree.fasta')

		 #infer lineage tree
                os.system('iqtree-omp -s ' + folder + cluster + '/tree.fasta -omp 3 -m TEST')

		#annotate with amino acid changes
                os.system('python /home/anna/Documents/software/for_anna/TreeAnnotateCmd.py ' + folder + cluster + '/seq_num.txt ' + folder + cluster + '/tree.fasta '\
                          + folder + cluster + '/tree.fasta.treefile ' + folder + cluster + '/CDR.txt ' + cluster + ' ' + folder + cluster + '/changes')
		#visualise tree
                tree_functions3.make_tree(folder + cluster + '/changes/' + cluster + '_annotated_treefile.new', folder + cluster + '/tree.png', clone_info)

     
                
if __name__=="__main__":
    main(sys.argv)

        

        
        
        
      
        
