'''

    Plot T-Aligner assembled mRNAs
    showing support for possitions

'''

import sys, argparse, matplotlib, os
from matplotlib import pyplot as plt
import matplotlib.patches as mpatch

from collections import defaultdict as ddict
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as mticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm



def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--r', '--fa',
        action='store',
        help="Reference fasta file (provide the same file that was used for mRNA assembly)")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output .pdf or .png file name")

    optionParser.add_argument('--t', '--taf',
        action='store',
        help="File with raw read alignments in TAF format (version 2)")

    optionParser.add_argument('--m', '--mrna',
        action='store',
        help="File in TAF format (version 2) with assembled mRNAs")

    optionParser.add_argument('--k', '--kmer',
        action='store',
        help="Split reads into T-less kmers for support calculation")

    optionParser.add_argument('--c', '--canonical',
        action='store',
        help="Number of canonical mRNA in TAF file, 1-based")

    return optionParser


def ReadFasta(input_lines, split_space = False):

    fasta_data = []

    for line in input_lines:
        line = line.rstrip()
        if len(line) < 1:
            continue
        if line[0] == '>':
            if split_space == False:
                fasta_data.append([line[1:],''])
            else:
                fasta_data.append([line[1:].split(' ')[0],''])
        else:
            fasta_data[-1][1] += line.lstrip().rstrip().upper()

    for i, entry in enumerate(fasta_data):
        fasta_data[i].append(len(entry[1]))

    return fasta_data


if __name__ == '__main__':

    runArgs = GetOptParser().parse_args(sys.argv[1:])

    inputFastaData = None
    if runArgs.r is not None:
        if os.path.isfile(runArgs.r):
            inputFastaData = open(runArgs.r, 'r').readlines()
    else:
        sys.stderr.write('No reference file found or -r option missing!')
        exit()

    fastaEntries = ReadFasta(inputFastaData)
    referenceSequence = fastaEntries[0][1]
    ref_tless_length = 0
    ref_tless_seq = ''

    for i in referenceSequence.lstrip().rstrip().upper():
        if i != 'T':
            ref_tless_length += 1
            ref_tless_seq = ref_tless_seq + i

    inputMrnaTaf = None
    if runArgs.m is not None:
        if os.path.isfile(runArgs.m):
            inputMrnaTaf = open(runArgs.m, 'r').readlines()
    else:
        sys.stderr.write('No mRNAs TAF file or -m option missing!')
        exit()


    plt.figure(figsize=[0.1 * ref_tless_length, 1 + 0.1 * len(inputMrnaTaf)])
    plt.gcf().set_dpi(300)
    ax = plt.gca()
    for i in range(ref_tless_length):
        nuc_rect = mpatch.Rectangle((10.0 + 10.0 * i, 0), 10.0, 10.0, edgecolor='#161717', facecolor="#34eba1", linewidth=1)
        ax.add_artist(nuc_rect)
        rx, ry = nuc_rect.get_xy()
        ax.annotate(ref_tless_seq[i], (rx+5.0, ry+5.0), color='black', fontsize=4, ha='center', va='center')

    inputReadsTaf = None
    if runArgs.t is not None:
        if os.path.isfile(runArgs.t):
            inputReadsTaf = open(runArgs.t, 'r').readlines()
    else:
        sys.stderr.write('No raw read mappings in TAF or -t option missing!')
        exit()

    canonicalMrnaNumber = -1
    if runArgs.c is not None:
        canonicalMrnaNumber = int(runArgs.c)

    kMerLength = -1
    kmerMode = False
    if runArgs.k is not None:
        kMerLength = int(runArgs.k)
        if (kMerLength % 2) != 0 or kMerLength > 11 or kMerLength < 3:
            sys.stderr.write('K-mer length must be even number in range (3-11)!')
            exit()
        else:
            kmerMode = True

    es_support_matrix = ddict(lambda : ddict(lambda : 1.0))
    es_tlp_cover = ddict(lambda : 1.0)

    # populate ES supports for simple mode
    if kmerMode == False:
        for taf_read_entry in inputReadsTaf:
            toks = taf_read_entry.rstrip().split('\t')
            p = int(toks[2])
            es = toks[4].split(';')[:-1]
            for i, es in enumerate(es):
                es_tlp_cover[p+i] += 1
                es_support_matrix[p+i][es] += 1
    # populate ES supports for Kmer mode
    else:
        for taf_read_entry in inputReadsTaf:
            toks = taf_read_entry.rstrip().split('\t')
            p = int(toks[2])
            es_vec = [str(es) for es in toks[4].split(';')[:-1]]
            for i in range(len(es_vec) - kMerLength + 1):
                es_tlp_cover[p+i+kMerLength/2] += 1
                es_support_matrix[';'.join(es_vec[i:i+kMerLength])][p+i+kMerLength/2] += 1

    color_norm = matplotlib.colors.Normalize(vmin=0, vmax=1.0, clip=True)
    color_mapper = cm.ScalarMappable(norm=color_norm, cmap=plt.get_cmap('viridis'))
    color_mapper_font = cm.ScalarMappable(norm=color_norm, cmap=plt.get_cmap('hot_r'))

    current_y  = 10.0
    currentTaf = 0
    max_supp_mrna_as = 0
    max_supp_mrna_ae = 0
    max_supp_mrna_value = 0.0

    for taf_mrna_entry in inputMrnaTaf:
        currentTaf +=1
        toks = taf_mrna_entry.rstrip().split('\t')
        #print('Mrna starts at {}, end at {}, first chars {}, ref chars {}'.format(toks[2], toks[3], toks[5][0:4], ref_tless_seq[int(toks[2]):int(toks[2])+4]))
        a_s, a_e = int(toks[2]), int(toks[3])
        e_pos = 0

        supp_mrna_value = 0.0
        if kmerMode == False:
            e_dr = toks[4].split(';')
            for i in range(a_s, a_e):
                rx = 15.0 + 10.0 * i; ry = current_y+5.0
                es = e_dr[e_pos]
                es_color = es_support_matrix[i][es] / es_tlp_cover[i]
                supp_mrna_value += es_color
                es_color_font = color_mapper_font.to_rgba(es_color)
                es_color = color_mapper.to_rgba(es_color)
                es_rect = mpatch.Circle((rx, ry), radius=5.0, edgecolor='#161717', facecolor=es_color, linewidth=0)
                ax.add_artist(es_rect)
                ax.annotate(es, (rx, ry), color=es_color_font, fontsize=4, ha='center', va='center')
                e_pos += 1

        else:
            e_dr = [str(es) for es in toks[4].split(';')[:-1]]
            e_pos += kMerLength/2
            for i in range(int(a_s + kMerLength/2), int(a_e - kMerLength/2) + 1):
                rx = 15.0 + 10.0 * i; ry = current_y+5.0
                es_kmer = ';'.join(e_dr[int(e_pos - kMerLength/2) : int(e_pos + kMerLength/2)])
                es_color = es_support_matrix[es_kmer][i] / es_tlp_cover[i]
                supp_mrna_value += es_color
                es_color_font = color_mapper_font.to_rgba(es_color)
                es_color = color_mapper.to_rgba(es_color)
                es_rect = mpatch.Circle((rx, ry), radius=5.0, edgecolor='#161717', facecolor=es_color, linewidth=0)
                ax.add_artist(es_rect)
                ax.annotate(e_dr[int(e_pos)], (rx, ry), color=es_color_font, fontsize=4, ha='center', va='center')
                e_pos += 1
        if max_supp_mrna_value < supp_mrna_value:
            max_supp_mrna_value = supp_mrna_value
            max_supp_mrna_y = current_y
            max_supp_mrna_as = a_s
            max_supp_mrna_ae = a_e

        if currentTaf == canonicalMrnaNumber:
            bounding_box_canonical = mpatch.Rectangle((8.0 + 10.0 * a_s, current_y), 10.0 * (a_e-a_s) + 2.0, 10.0, facecolor=(1,1,1,0), edgecolor='#ff1a1a', linewidth=1.0)
            ax.add_artist(bounding_box_canonical)
        current_y += 10.0

    bounding_box_realbest = mpatch.Rectangle((8.0 + 10.0 * max_supp_mrna_as, max_supp_mrna_y), 10.0 * (max_supp_mrna_ae-max_supp_mrna_as) + 2.0, 10.0, facecolor=(1,1,1,0), edgecolor='#a434eb', linewidth=2.0)
    ax.add_artist(bounding_box_realbest)

    ax.set_xlim((0, ref_tless_length * 10.0 + 20.0))
    ax.set_ylim((0, len(inputMrnaTaf) * 10.0 + 10.0))
#
#   END PLOTTING COMPASS
#
    output_pic = 'output.pdf'
    if runArgs.o is not None:
        output_pic = str(runArgs.o)
    plt.savefig(output_pic)

#python mrnaplotter.py --r RPS12__SylvioX10.consensus.fasta --m RPS12_25_assembled_mrna.taf --t SYL_taf/SYL_RPS12__SylMito__REP1_mapped_reads.taf --o output.png
