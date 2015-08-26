#miscallaneous helper functions to be used by scripts which analyse CDR3s

import re

def create_list(infile):
"""creates a list of (folder,fraction) tuples from a text file with tab separated folder_names and fractions"""
    with open(infile) as file_list:
        files = file_list.readlines()
        files2 = [line.strip().split() for line in files]
        return files2


def all_seqs(sequence):
"""can be used where a condition function is required but it is desired to analyse all sequences"""
    return True


def no_cysteines(sequence):
"""returns TRUE if an amino acid sequence contains no cysteines"""
    if sequence.count('C') == 0:
        return True
    else:
        return False


def gt1_cysteines(sequence):
"""returns TRUE if an amino acid sequence contains two or more cysteines"""
    try:
        if sequence.count('C') > 1:
            return True
        else:
            return False
    except:
        return False


def non_canonical(sequence):
"""returns TRUE if an amino acid sequence contains two or more cysteines with a 
non-canonical spacing (i.e. not 3 or 4)"""
    if sequence.count('C') >1:
        pat = re.compile('C....?C')
        info = pat.search(sequence)
        if info:
            return False
        else:
            return True
    return False
    

def spacing_3(sequence):
"""returns TRUE if an amino acid sequence contains two or more cysteines with a 
spacing of 3 amino acids"""
    if sequence.count('C') >1:
        pat = re.compile('C...C')
        info = pat.search(sequence)
        if info:
            return True
    return False


def spacing_4(sequence):
"""returns TRUE if an amino acid sequence contains two or more cysteines with a 
spacing of 4 amino acids"""
    if sequence.count('C') >1:
        pat = re.compile('C....C')
        info = pat.search(sequence)
        if info:
            return True
    return False


def find_v_gene(v_gene):
"""returns V gene if IMGT all assignments are in agreement; allele information is not considered"""
    v_genes = set()
    pat = re.compile('IGHV.')
    for item in pat.finditer(v_gene):
        v_genes.add(item.group(0))
    if len(v_genes) == 1:
        return v_genes.pop()
    else:
        return False
            

def nc_macaque(sequence):
"""returns TRUE if an amino acid sequence contains two or more cysteines with a 
non-canonical spacing for macaques (i.e. 4)"""
    if sequence.count('C') >1:
        pat = re.compile('C....C')
        info = pat.search(sequence)
        if info:
            return False
        else:
            return True
    return False   
        


