#functions to simulate the occurence of two or more cysteines in CDR3s

import random

def generate_random_seqs(lengths,num):
"""this simulation assumes all amino acids are equally likey to occur"""
    amino_acids = ['P','G','A','V', 'L', 'I', 'F', 'Y', 'W',\
            'M', 'C', 'S', 'T', 'N', 'Q', 'H', 'K', 'R', 'D', 'E']
    dat = []
    n = 0
    for i in lengths:
        dat.append([])
        for j in range(num):
            cdr = []
            for k in range(i):
                cdr.append(amino_acids[random.randrange(20)])
            dat[n].append(''.join(cdr))
        n += 1
    dat1 = []
    num2 = 0
    for item in dat:
        n = 0
        for seq in item:
            if seq.count('C')> 1:
                n += 1
        dat1.append([lengths[num2],n/float(num)])
        num2 += 1
    return dat1

def generate_random_seqs_2(lengths,num):
"""this simulation uses the codon based probability (3.28%) of cysteine occurence"""
    dat = []
    n = 0
    for i in lengths:
        dat.append([])
        for j in range(num):
            cdr = []
            for k in range(i):
                r_num = random.randrange(10000)
                if r_num < 328:
                    amino_acid = 'C'
                else:
                    amino_acid = 'X'
                cdr.append(amino_acid)
            dat[n].append(''.join(cdr))
        n += 1
    dat1 = []
    num2 = 0
    for item in dat:
        n = 0
        for seq in item:
            if seq.count('C')> 1:
                n += 1
        dat1.append([lengths[num2],n/float(num)])
        num2 += 1
    return dat1

def generate_random_seqs_3(lengths,num):
"""this simulation uses the frequency cysteines occur in mammalian proteins
to model cysteine occurence"""
    dat = []
    n = 0
    for i in lengths:
        dat.append([])
        for j in range(num):
            cdr = []
            for k in range(i):
                r_num = random.randrange(10000)
                if r_num < 226:
                    amino_acid = 'C'
                else:
                    amino_acid = 'X'
                cdr.append(amino_acid)
            dat[n].append(''.join(cdr))
        n += 1
    dat1 = []
    num2 = 0
    for item in dat:
        n = 0
        for seq in item:
            if seq.count('C')> 1:
                n += 1
        dat1.append([lengths[num2],n/float(num)])
        num2 += 1
    return dat1

