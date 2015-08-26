#miscellaneous functions for use in tree_pipeline

import re
import math
from ete2 import Tree, TreeStyle, TextFace,  NodeStyle, faces, AttrFace, CircleFace

def write_fasta_files(cluster, clusters, seqs, outfile1, outfile2):
    with open(outfile1, 'w') as first_timepoint:
        with open(outfile2, 'w') as all_seqs:
            pat = re.compile('(\*)+Cluster\s' + cluster + '(\*)+\n(.+?)(\*|\Z)', re.S)
            info = pat.search(clusters)
            first = True
            iden = False
            time_1 = ''
            time = ''
            text = info.group(0).split('\n')
            for line in text:
                if line and line[0] == '>':
                    if first:
                        time_1 = line.strip().split(';')[-1]
                        first = False
                    pat1 = re.compile('(' + line.strip() + '.*?\n\w+?)(\n|\Z)>?')
                    seq = pat1.search(seqs)
                    if line.strip().split(';')[-1] == time_1:
                        first_timepoint.write(seq.group(1) + '\n')
                    all_seqs.write(seq.group(1) + '\n')
                    

def process_fasta_for_tree(infile, outfile):
    with open(infile) as aligned:
        with open(outfile, 'w') as out:
            clone_info = {}
            out.write(aligned.readline())
            out.write(aligned.readline())
            for line in aligned:
                if line[0] == '>':
                    pat = re.compile('>clone-(.*?);(.*?);(.*?);(.+?)')
                    info = pat.match(line)
                    out.write('>' + info.group(1) + '  \n')
                    if info.group(1) in clone_info.keys():
                        clone_info[info.group(1) + '-' + info.group(4).strip()] = (info.group(4).strip(),info.group(2),info.group(3))
                    else:
                        clone_info[info.group(1)] = (info.group(4),info.group(2),info.group(3))
                else:
                    out.write(line)
            for key in clone_info.keys():
                print key, '\t', clone_info[key], '\n'
    return(clone_info)


def style_node(node, colour, size):
    nstyle = NodeStyle()
    nstyle["fgcolor"] = colour
    nstyle["size"] = size
    node.set_style(nstyle)


def style_node2(node, colour, size):
    C = CircleFace(radius = size, color = colour, style="sphere")
    C.opacity = 0.8
    node.add_face(C, 0, position="float")

    
def make_tree(treefile, image_file, clone_info):
    colour_list = ['MidnightBlue','RoyalBlue', 'LightSkyBlue', 'Aquamarine', 'SpringGreen', 'GreenYellow',\
                   'Gold','DarkOrange']
    weeks = ['16', '30', '38', '48', '59', '119', '176', '206']
    weeks = ['6', '14', '53', '92','144']
    t = Tree(treefile,format = 1)
    ts = TreeStyle()
    for i in range(5):
        ts.legend.add_face(CircleFace(20, colour_list[i]), column=0)
        ts.legend.add_face(TextFace('week' + weeks[i]), column=1)
    ts.legend_position = 2
    ts.show_leaf_name = True
    ts.branch_vertical_margin = 15
    ts.rotation = 90
    ns = NodeStyle()
    ns["size"] = 1
    ns.hz_line_width = 10
    ns.vt_line_width = 10
    edge = 0
    for node in t.traverse():
        node.name = node.name.replace("'", "")
        node.name = node.name.replace(".", ",")
        name = node.name.split(' ')[0]
        print name
        if name in clone_info.keys():
            style_node(node, colour_list[int(clone_info[name][0])-1], int(int(clone_info[name][1])/10)+5)
        if not node.is_leaf() and node.name != 'NoName':
                f = TextFace(node.name)
                f.margin_top = 2.5
                f.margin_bottom = 2.5
                f.margin_right = 2.5
                f.margin_left = 2.5
                node.add_face(f, column=0, position="branch-top")
    t.render(image_file, tree_style = ts)

    
            
    
