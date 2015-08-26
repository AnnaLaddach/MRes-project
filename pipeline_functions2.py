#miscellaneous functions for use in pipelines

import re
import glob
import string


def vidjil_num_seq(infile1, infile2, outfile1, outfile2):
    with open(infile1) as clone_file:
        clones = clone_file.read()
    with open(infile2) as seq_file:
        seqs = seq_file.read()
    out1 = open(outfile1, 'w')
    out2 = open(outfile2, 'w')
    pat1 = re.compile('>(clone-\d+).*?window\s*?\n([A-Z]+)')
    for match in pat1.finditer(clones):
        iden = match.group(1)
        window = match.group(2)
        num_seq = 0
        pat2 = re.compile('>\d+-(\d+)\s*?\n[A-Z]*?'+ window + '[A-Z]*?' )
        for match1 in pat2.finditer(seqs):
            num_seq += int(match1.group(1))
        window2 = reverse_complement(window)
        pat3 = re.compile('>\d+-(\d+)\s*?\n[A-Z]*?'+ window2 + '[A-Z]*?')
        for match1 in pat3.finditer(seqs):
            num_seq += int(match1.group(1))
        out1.write(iden + ' ' + str(num_seq) + '\n')
        out2.write('>' + iden + '\n' + window + '\n')


def reverse_complement(seq):
    complement_table = string.maketrans('ATCG','TAGC')
    seq1 = seq[::-1]
    return(seq1.translate(complement_table))


def remove_single_reads(infile, outfile):
    with open(infile) as file:
        out = open(outfile, 'w')
        for line in file:
            if line[0] == '>':
                if line.strip()[-2::] !='-1':
                    not_single = True
                    out.write(line)
                else:
                    not_single = False
            else:
                if not_single:
                   out.write(line)


def extract_CDR3(infile, outfile):
    with open(infile) as file:
        out = open(outfile, 'w')
        file.readline()
        for line in file:
            attributes = line.split('\t')
            if attributes[2] == 'productive' and attributes[2] != 'unproductive':
                try: 
                    out.write('>' + attributes[1] + '\n' + attributes[14] + '\n')
                except IndexError:
                    print 'no results for ' + attributes[1]


def produce_VDJ_summary_vidjil(infile, outfile):
    with open(infile) as file:
        out = open(outfile, 'w')
        seqs = file.read()
        pat = re.compile('>(clone.\d+).+?(IGHV[^\s]+)\s.+?(IGHD[^\s]+)\s.+?(IGHJ[^\s]+)') 
        for match in pat.finditer(seqs):
            out.write(match.group(1) + ',' + match.group(2) + ',' + match.group(3) +',' + match.group(4) + ',\n')


def produce_VDJ_summary_imgt(summary, n_sequence, a_sequence, outfile):
    with open(summary) as summary_file:
        with open(n_sequence) as n_file:
            with open(a_sequence) as a_file:
                summary_file.readline()
                n_file.readline()
                a_file.readline()
                out = open(outfile, 'w')
                for s_line in summary_file:
                    attributes = s_line.split('\t')
                    nucs = n_file.readline().split('\t')
                    aas = a_file.readline().split('\t')
                    clone = find_clone_name(attributes[1])
                    try:
                        v_seg = find_segs('IGHV', attributes[3])
                        d_seg = find_segs('IGHD', attributes[13])
                        j_seg = find_segs('IGHJ', attributes[9])
                        out.write(clone + ',' + '-or-'.join(v_seg) + ',' + '-or-'.join(d_seg) + ',' + '-or-'.join(j_seg)\
                                    + ',' + nucs[14] + ',' + aas[14] + ',' + nucs[6] + ',' + aas[9] + ',' + aas[10] + ',' +\
                                    aas[11] + ',' + aas[12] + ',' + aas[13] + ',' + aas[17] + ',' + attributes[14] +',\n')
                    except IndexError:
                        print 'no results for ' + clone
                        
                    
def find_segs(seg_type, line):
    pat = re.compile('(' + seg_type + '.*?)(\*|[A-Z])')
    segs = set()
    try:
        for match in pat.finditer(line):
            seg = match.group(1)
            if len(seg) > 8:
                seg = seg[:8]
            segs.add(seg)
    except AttributeError:
        segs = 'none'
    return segs
      

def determine_prevailing_germline_seg(clusters, germline, numbers, outfile):
    with open(numbers) as num_file:
        num = num_file.read()
    with open(germline) as germ_file:
        germ = germ_file.read()
    out = open(outfile, 'w')
    first = True
    cluster_num = 0
    with open(clusters) as cluster_file:
        for line in cluster_file:
            if line[0] == '>':
                if not first:
                    find_representative_gene_segs(num_germline_segments, cluster_num, out)
                    cluster_num += 1
                num_germline_segments = dict()
                first = False
            elif line[0] != '>':
                clone = find_clone_name(line)
                gene_segs = search_for_germline_genesegments(clone, germ)
                number = find_frequency_of_clone(clone, num)
                if gene_segs in num_germline_segments.keys():
                    num_germline_segments[gene_segs] += int(number)
                else:
                    num_germline_segments[gene_segs] = int(number)
        find_representative_gene_segs(num_germline_segments, cluster_num, out)


def find_clone_name(line):
    pat = re.compile('(clone.\d+)')
    info = pat.search(line)
    return info.group(1)


def search_for_germline_genesegments(clone, germ):
    pat = re.compile(clone + '.*?,(.*?,.*?,.*?),.*?,.*?,.*?,.*?,.*?,.*?,.*?')
    info = pat.search(germ)
    try:
        gene_segs = info.group(1)
    except AttributeError:
        gene_segs = 'none'
        print 'no gene segments found for ' + clone
    return gene_segs


def search_for_segs_info(clone, germ):
    pat = re.compile(clone + '.*?,(.*?,.*?,.*?),(.*?,.*?,.*?,.*?,.*?,.*?,.*?,.*?,.*?,.*?),')
    info = pat.search(germ)
    segs = info.group(1)
    cdr = info.group(2)
    return(segs, cdr)


def find_frequency_of_clone(clone, num):
    pat = re.compile(clone + '\s(\d+)\s')
    info = pat.search(num)
    number = info.group(1)
    return number


def find_representative_gene_segs(num_germline_segments, cluster_num, out):
    try:
        representative = max(num_germline_segments, key = num_germline_segments.get)
        out.write(representative + ',' + str(num_germline_segments[representative]) + '\n')
    except ValueError:
        print 'no prevailing germline gene-segments found for cluster' + str(cluster_num)
        out.write('\n')


def find_germline_prevail(p_line):
    pat = re.compile('(.*?,.*?,.*?),(\d+)')
    info = pat.match(p_line)
    try:
        segs1 = info.group(1)
        number = info.group(2)
    except AttributeError:
        print 'no prevailing germline gene segments found for cluster'
        segs1 = 'none'
        number = -1
    return (segs1, number)


def find_sequence(clone, seqs):
    pat = re.compile('>' + clone + '-.*?\n([A-Z]+.*?)(>|\Z)', re.S)
    seq_info = pat.search(seqs)
    spaced_seq = seq_info.group(1)
    space = re.compile('\s')
    nucleotides = space.sub('', spaced_seq)
    return nucleotides

    
def summary_w(summary_file, total, summary, IMGTflag):
    if IMGTflag:
        summary_file.write('cluster,v-segment,d-segment,j-segment,frequency,relative frequency,CDR3 nucleotides,CDR3 amino acids,VDJ nucleotides,FR1 amino acids,'\
                           + 'CDR1 amino acids,FR2 amino acids,CDR2 amino acids,FR3 amino acids,FR4 amino acids,D reading frame,representative sequence,total\n')
        for item in summary:
            summary_file.write(item[0] + ',' + item[1] + ',' + str(item[2]) + ',' + str(float(item[2])/total) + ',' + item[3] + ',' + item[4] + ',' + str(total) + '\n')
    else:
        summary_file.write('cluster, v-segment, d-segment, j-segment, frequency, relative frequency, representative sequence,\n')
        for item in summary:
            summary_file.write(item[0] + ',' + item[1] + ',' + str(item[2]) + ',' + str(float(item[2])/total) + ',' + item[3] + '\n')


def create_summary_elements(cluster_num, segs1, number, nucleotides, cdr_info, IMGTflag):
    if IMGTflag:
        return ['Cluster ' + str(cluster_num), segs1, number, cdr_info, nucleotides]
    else:
        return ['Cluster ' + str(cluster_num), segs1, number, nucleotides]
    

def write_summary(clusters, prevailing, germline, seq, numbers, summary_f, correct_clusters, IMGTflag = True):
    with open(seq) as sequence_file:
        seqs = sequence_file.read()
    with open(numbers) as num_file:
        num_clones = num_file.read()
    with open(germline) as germ_file:
        germ = germ_file.read()
    summary = []
    extended_summary = []
    first = True
    summary_file = open(summary_f, 'w')
    cluster_file = open(correct_clusters, 'w')
    cluster_num = 0
    with open(prevailing) as prevail_file:
        with open(clusters) as clust:
            for line in clust:
                if line[0] == '>':
                    if not first:
                        nucleotides = find_sequence(rep_clone, seqs)
                        summary_elements = create_summary_elements(cluster_num, segs1, number, nucleotides,\
                                                                   cdr_info, IMGTflag)
                        summary.append(summary_elements)
                        cluster_num += 1 
                    rep_clone = ''
                    cdr_info = ''
                    num = 0
                    p_line = prevail_file.readline()
                    (segs1, number) = find_germline_prevail(p_line)
                    cluster_file.write(line)
                    first = False
                elif line[0] != '>':
                    clone = find_clone_name(line)
                    if IMGTflag:
                        (segs2, extra_info) = search_for_segs_info(clone, germ)
                    else:
                        segs2 = search_for_germline_genesegments(clone, germ)
                    if segs1 == segs2:
                        cluster_file.write(line)
                        num1 = find_frequency_of_clone(clone, num_clones)
                        if int(num1) > num:
                            num = int(num1)
                            rep_clone = clone
                            if IMGTflag:
                                cdr_info = extra_info
            nucleotides = find_sequence(rep_clone, seqs)
            summary_elements = create_summary_elements(cluster_num, segs1, number, nucleotides,\
                                                                   cdr_info, IMGTflag)
            summary.append(summary_elements)
            summary = sorted(summary, key = lambda frequency: int(frequency[2]), reverse = True)
            total = sum([int(item[2]) for item in summary])
            summary_w(summary_file, total, summary,IMGTflag)
            return total


def determine_cluster(line):
    pat = re.compile('(Cluster\s\d+),')
    info = pat.search(line)
    return info.group(1)


def find_clones(cluster, clean):
    pat = re.compile('>'+ cluster + '.+?(>Cluster|\Z)', re.S)
    info = pat.search(clean)
    return info.group(0)


def write_sequence_summary(summary_file, clean_clusters, seq, numbers, out, total, IMGTflag = True, germline = ''):
    out_file = open(out, 'w')
    out_file2 = open(out + '_g', 'w')
    out_file2.write('clone,v-segment,d-segment,j-segment,frequency,relative frequency,CDR3 nucleotides,CDR3 amino acids,VDJ nucleotides,FR1 amino acids,'\
                           + 'CDR1 amino acids,FR2 amino acids,CDR2 amino acids,FR3 amino acids,FR4 amino acids,D reading frame,sequence\n')
    if IMGTflag:
        with open(germline) as germ_file:
            germ = germ_file.read()
    with open(numbers) as num_file:
        num_clones = num_file.read()
    with open(seq) as sequence_file:
        seqs = sequence_file.read()
    with open(clean_clusters) as clusters:
        clean = clusters.read()
    with open(summary_file) as summary:
        summary.readline()
        for line in summary:
            cluster = determine_cluster(line)
            out_file.write('>' + cluster + '\n')
            summary_info = line.split(',')
            clones = find_clones(cluster, clean)
            pat = re.compile('clone.\d+')
            clone_info = []
            for match in pat.finditer(clones):
                clone = match.group(0)
                nucleotides = find_sequence(clone, seqs)
                number = find_frequency_of_clone(clone, num_clones)
                if IMGTflag:
                    CDR_info = search_for_segs_info(clone, germ)[1]
                    clone_info.append([clone, int(number), float(number)/total, CDR_info, nucleotides])
                else:
                    clone_info.append([clone, int(number), float(number)/total, nucleotides])
            clone_info = sorted(clone_info, key = lambda frequency: frequency[2], reverse = True)
            if IMGTflag:
                for item in clone_info:
                    out_file.write(item[0] + ',' + str(item[1]) + ',' + str(item[2]) + ',' + item[3] + ',' + item[4] + '\n')
                    out_file2.write(item[0] + ',' + summary_info[1] + ',' + summary_info[2] + ',' + summary_info[3] + ',' + str(item[1]) + ',' + str(item[2]) + ',' + item[3] + ',' + item[4] + '\n')
            else:
                for item in clone_info:
                    out_file.write(item[0] + ',' + str(item[1]) + ',' + str(item[2]) + ',' + item[3] + '\n')
    
             
            
def split_vidjil_output(infile, outfile):
    with open(infile) as file:
        f = 1
        n = 1
        out = open(outfile + str(f) + '.fasta', 'w')
        for line in file:
            if line[0] == '>':
                n +=1
            if n % 100000 == 0:
                out = open(outfile + str(f) + '.fasta', 'w')
                f += 1
            out.write(line)
        return f

                
def combine_IMGT_output(foldername, filename, num_files, outfile):
    out = open(outfile, 'w')
    for i in range(1, num_files + 1):
        first = True
        with open(glob.glob(foldername + str(i) + '/' + filename + '*')[0]) as file:
            if i == 1:
                out.write(file.read())
            else:
                n = 1
                for line in file:
                    if n != 1:
                        out.write(line)
                        n += 1


def write_all_seq_summary(infile, clone_num, outfile, total):
    with open(clone_num) as num:
        numbers = clone_num.read()
    with open(outfile, 'w') as out:
        with open(infile) as nucs:
            for line in nucs:
                attributes = line.split('\t')
                if attributes[2] == 'productive' and attributes[2] != 'unproductive':
                    pass
                    
               




    
            
                            
                        
   
                    
                    
    
        











                    
                
