'''

    Plot T-Aligner output file in TAF2 format
    as compass-plot

'''

import sys, argparse, matplotlib, os
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as mticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm



def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--r', '--fa',
        action='store',
        help="Reference fasta file (first reference sequence will be used)")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output pdf/png file name")

    optionParser.add_argument('--t', '--taf',
        action='store',
        help="File with read alignments in TAF format (version 2)")

    optionParser.add_argument('--c', '--cthr',
        action='store',
        help="Editing state coverage threshold [int, default=10]")


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

    for i in referenceSequence.lstrip().rstrip().upper():
        if i != 'T':
            ref_tless_length += 1

    compass_matrix = []

    for y in range(-20, 20):
        compass_matrix.append([[0 for x in range(-20, 20)] for x in range(ref_tless_length)])

    # read TAF

    inputTafFile = ''
    if runArgs.t is not None:
        if os.path.isfile(runArgs.t):
            inputTafFile = runArgs.t
    else:
        sys.stderr.write('No TAF alignments file found or -t option missing!')
        exit()

    with open(inputTafFile, 'r') as ifile:
        for line in ifile:
            toks = line.split('\t')
            s = int(toks[2])
            es = [int(x) for x in toks[4][:-1].split(';')[1:-2]]

            # skip VERY long insertions
            skip = False
            for ess in es:
                if ess > 20:
                    skip = True
                    break
            if skip:
                continue
            # populate compass matrix
            for p, ess in enumerate(es):
                if s + p >= ref_tless_length or p == 0:
                    continue
                compass_matrix[20-ess][s + p][20-es[p-1]] += 1

    es_cover_thr = 10
    if runArgs.c is not None:
        es_cover_thr = int(runArgs.c)

#
#   PLOTTING COMPASS
#
    plt.figure(figsize=[ref_tless_length / 6, 12])
    plt.gcf().set_dpi(300)
    plt.plot([-10.0, 10.0 * ref_tless_length + 10.0],[210,210], c='#a1d76a', lw=10, alpha=0.3)

    ax = plt.gca()
    x_ticks_list = [x for x in range(10 * 25, 10 * ref_tless_length - 100, 10 * 25)]
    x_ticks_labels = ['{} nt'.format(int(x/10.0)) for x in x_ticks_list]

    ax.set_xticks(x_ticks_list)
    ax.set_xticklabels(x_ticks_labels)
    plt.xticks(fontsize=36, fontweight='bold')

    ax.grid(linestyle='-', linewidth=1, which='both')

    plt.xlim(-10.0, 10.0 * ref_tless_length + 10.0)
    ax.set_yticks((90, 120, 150, 180, 210, 240, 270, 300, 330))
    ax.set_yticklabels(('-12 U', '-9 U', '-6 U', '-3 U', 'ref', '+3 U', '+6 U', '+9 U', '+12 U'))
    plt.yticks(fontsize=36, fontweight='bold')
    plt.ylim(110, 340)
    
    for y in range(len(compass_matrix)):
        for x in range(len(compass_matrix[0])):

            vectors = [[n, m] for n, m in zip(list(range(len(compass_matrix[y][x]))), compass_matrix[y][x])]
            vectors.sort(key = lambda k : k[1])
            xc = 10.0 * x + 5.0
            yc = 410 - 10.0 * y - 5.0

            # look at last 3 vector components
            ins_p = [-9.0, -6.0, -4.0]
            del_p = [ 9.0,  6.0,  4.0]
            color_p = ['blue', 'red', 'black']
            for i in range(-1, -4, -1):
                if vectors[i][1] < es_cover_thr or (vectors[-1][1] - vectors[i][1]) > 0.2 * vectors[-1][1]:
                    continue
                dr = vectors[i][0]
                es = 20 - dr

                y_shift = 0
                if es < 0:
                    y_shift = ins_p.pop()
                elif es > 0:
                    y_shift = del_p.pop()
                yi = 410 - 10.0 * dr - 5.0
                plt.arrow(xc, yc, -6.0, y_shift, width=0.1, color=color_p.pop())




#
#   END PLOTTING COMPASS
#
    output_pic = 'output.pdf'
    if runArgs.o is not None:
        output_pic = str(runArgs.o)
    plt.savefig(output_pic)
