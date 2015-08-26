#Script to cluster sequences(i.e. CDR3s) based on fraction identity after Needleman-Wunsch alignment using single-linkage clustering.
#Usage: python cluster.py seq.fasta threshold_identity max_aa_indels indel_penalty
#seq.fasta - input sequences file
#threshold_identity - minimum proportion identity by which to cluster (from 0 to 1)
#max_aa_indels - maximum number of amino acid indels from one link to the next to accept when clustering (although any number can be specified it is assumed a maximum of 0-2 will be relevant in a biological context). 
#indel_penalty - number to be subtracted from proportion identity if indels are present e.g. 0.1 (this does not increase if more than one insertion or deletion is allowed).
#An output file will be created with the name seq_threshold_identity_max_aa_indels_indel_penalty.txt, where seq is the input filename without the fasta extension.  

from Bio import pairwise2
import sys

class Seq:
    def __init__(self, iden, seq):
        """creates a sequence object with a sequence (.seq) and identifier (.iden)"""
        self.seq = seq
        self.iden = iden
        self.length = len(seq)


class Cluster:
    def __init__(self, iden, seqs):
        """creates a cluster object with a list of sequences (.seq) and an identifier (.iden)"""
        self.seqs = [seqs]
        self.iden = iden

    def add_seq(self, seq):
        """adds a sequence to the cluster's sequence list"""
        self.seqs.append(seq)


def read_fasta(in_file):
    """reads a fasta file and returns a list of sequence objects sorted by sequence length"""
    with open(in_file) as file:
        nucs = ''
        iden = ''
        sequences = []
        for line in file:
            if line[0]== '>':
                if nucs:
                    sequences.append(Seq(iden, nucs))
                iden = line[1:].strip()
                nucs = ''
            else:
                nucs += line.strip()
        sequences.append(Seq(iden, nucs))
        sequences.sort(key = lambda x: x.length, reverse = True)
        return sequences


def create_clusters(sequences, threshold, max_aa_indels, indel_penalty):
    """clusters sequences using single linkage and specified parameters given a list of sequence objects.\
Returns a list of cluster objects"""
    clusters = []
    n = 0
    num_added = 0
    while sequences: #iterate until all sequences have been assigned a cluster.
        print 'cluster', n
        length = sequences[0].length
        clusters.append(Cluster(n,sequences[0])) #initiate a new cluster.
        sequences = sequences[1:] #remove sequence assigned to new cluster from sequence list.
        first = True
        while num_added > 0 or first: #if no new sequences have been added during the previous iteration the cluster has been completed.
            num_added = 0
            seqs_to_add = [] #initiate list of sequences to be added to cluster.
            for sequence in list(sequences):#double iteration through sequences in a copy of sequence list (to prevent problems with deleting list items) and sequences in current cluster.
                if length - sequence.length > max_aa_indels:
                    break
                for seq in clusters[n].seqs:
                    if align(sequence, seq, max_aa_indels, indel_penalty) > threshold:
                        seqs_to_add.append(sequence) #if sequence from sequence list passes threshold add it to seqs_to_add list and remove it from the sequence list.
                        length = sequence.length
                        num_added += 1
                        sequences.remove(sequence)
                        break #move on to the next sequence to prevent a sequence being assigned to the same cluster multiple times.
            if seqs_to_add:
                for seq in seqs_to_add: #add sequences_to_add to the sequence list of cluster.
                    clusters[n].add_seq(seq)
            first = False
        n += 1
    return clusters


def align(seq1, seq2, max_aa_indels, indel_penalty):
    """aligns sequences of two sequences using the Needleman-Wunsch algorithm, returns the fraction\
identitity of the shorter sequence unless indels exceed the max_aa_indels parameter specified. In this case 0 is returned."""
    if abs(len(seq1.seq) - len(seq2.seq)) > max_aa_indels*3:
           return 0
    elif len(seq1.seq) > len(seq2.seq):
           min_len = len(seq2.seq)
    else:
           min_len = len(seq1.seq) 
    alignment = pairwise2.align.globalxx(seq1.seq, seq2.seq)
    n = 0
    for i in range(len(alignment[0][0])):
        if alignment[0][0][i] != '-' and alignment[0][0][i] == alignment[0][1][i]:
            n += 1
    identity = float(n)/min_len
    if len(seq1.seq) != len(seq2.seq):
        identity -= indel_penalty #if there is an indel decrease the identity by the user specified amount.
    return identity 


def main(argv):
    """Reads and clusters sequences in a fasta file using single linkage and user defined parameters"""
    if len(argv) != 5 :
         print 'usage: python cluster.py seq.fasta threshold_identity max_aa_indels indel_penalty'
         sys.exit()
    sequences = read_fasta(argv[1])
    clusters = create_clusters(sequences, float(argv[2]), int(argv[3]), float(argv[4]))
    with open(argv[1][:-5] + '_cluster_' + argv[2] + '_' + argv[3] + '_' + argv[4] + '.txt', 'w') as out_file:
        for cluster in clusters:
            out_file.write('******Cluster ' + str(cluster.iden) + '******\n')
            for seq in cluster.seqs:
                out_file.write('>' + seq.iden + '\n')
                out_file.write(seq.seq + '\n')
            out_file.write('\n')


if __name__=="__main__":
    main(sys.argv)

    

    
                    
    
        
        



            
