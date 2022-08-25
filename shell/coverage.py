import sys, argparse, matplotlib, os
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as mticker

# INPUT : - reference sequence (cryptogene)
#         - alignments in TAF format
#         - figure size (reference start:end)
def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--r', '--fa',
        action='store',
        help="Reference fasta file [needs reference name if multifasta is given]")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output pdf/png file name")

    optionParser.add_argument('--l', '--loc',
        action='store',
        help="Locus of reference to plot [1-based reference coordinates]")

    optionParser.add_argument('--t', '--taf',
        action='store',
        help="File with read alignments in TAF format")

    optionParser.add_argument('--n', '--ref',
        action='store',
        help="Name of the fasta record")

    optionParser.add_argument('--g', '--gt',
        action='store',
        help="Type of the graph: bars, dots, prod")


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
    inputFastaName = 0
    if runArgs.n is not None:
        for i, n in enumerate(fastaEntries):
            if n[0] == runArgs.n:
                inputFastaName = i
                break

    referenceSequence = fastaEntries[inputFastaName][1]
    # get cryptogene reference length
    ref_tless_length = 0

    for i in referenceSequence.lstrip().rstrip().upper():
        if i != 'T':
            ref_tless_length += 1

    # fill matrices with zero values
    es_read_matrix = []
    es_productive_matrix = []

    for y in range(-20, 20):
        es_read_matrix.append([0 for x in range(ref_tless_length)])
        es_productive_matrix.append([-0.01 for x in range(ref_tless_length)])

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

            # detect if editing drops to reference (cutoff 7)
            for p, ess in enumerate(es):
                if s + p >= ref_tless_length:
                    continue
                es_read_matrix[20-ess][s + p] += 1
                if p > 7:
                    if es[p-8:p] == [0, 0, 0, 0, 0, 0, 0, 0]:
                        es_productive_matrix[20-ess][s + p] += 1

    l_from = 0
    l_to   = ref_tless_length
    if runArgs.l is not None:
        f, t = str(runArgs.l).split(',')
        l_from = int(f)
        l_to   = int(t)

    #plt.figure(figsize=[(l_to-l_from) / 8, 6])
    plt.figure(figsize=[ref_tless_length / 6, 12])
    plt.gcf().set_dpi(300)

    graph_type = 'bars'
    if runArgs.g is not None:
        graph_type = str(runArgs.g)

    if graph_type == 'prod':
        xi = []
        yi = []
        zi = []
        for y in range(len(es_productive_matrix)):
            for x in range(len(es_productive_matrix[0][l_from:l_to])):
                xi.append(10.0 * x)
                yi.append(410 - 10.0 * y)
                if y != 20 and es_read_matrix[y][x] > 20:
                    zi.append(es_productive_matrix[y][x] / (1.0+es_read_matrix[y][x]))
                else:
                    zi.append(0)

        plt.plot([10.0 * l_from - 10.0, 10.0 * l_to + 10.0],[210,210], c='#a1d76a', lw=10, alpha=0.3)
        plt.scatter(xi,yi, c=zi, cmap='Oranges', s=10.0)
        #plt.plot(main_orf_X, main_orf_Y, color='#a1d76a', lw=3)


    if graph_type == 'dots':
        xi = []
        yi = []
        zi = []
        for y in range(len(es_read_matrix)):
            for x in range(len(es_read_matrix[0][l_from:l_to])):
                xi.append(10.0 * x)
                yi.append(410 - 10.0 * y)
                if y != 20 and es_read_matrix[y][x] > 20:
                    zi.append(es_read_matrix[y][x])
                else:
                    zi.append(0)

        plt.plot([10.0 * l_from - 10.0, 10.0 * l_to + 10.0],[210,210], c='#a1d76a', lw=10, alpha=0.3)
        scatter_tal = plt.scatter(xi,yi, c=zi, cmap='Greens', s=80.0, norm=matplotlib.colors.LogNorm()).set_zorder(20)

        ax = plt.gca()
        x_ticks_list = [x for x in range(10 * 25, 10 * ref_tless_length - 100, 10 * 25)]
        #x_ticks_labels = ['{}'.format(int(x/10.0)) for x in x_ticks_list]; x_ticks_labels[-1] = x_ticks_labels[-1] + ' nonT'
        x_ticks_labels = ['{} nt'.format(int(x/10.0)) for x in x_ticks_list]

        ax.set_xticks(x_ticks_list)
        ax.set_xticklabels(x_ticks_labels)
        plt.xticks(fontsize=36, fontweight='bold')

        ax.grid(linestyle='-', linewidth=3, which='both')

        plt.xlim(-10.0, 10.0 * ref_tless_length + 10.0)
        #ax.set_ylabel('Editing state',fontsize=25, fontweight='bold')
        ax.set_yticks((90, 120, 150, 180, 210, 240, 270, 300, 330))
        ax.set_yticklabels(('-12 U', '-9 U', '-6 U', '-3 U', 'ref', '+3 U', '+6 U', '+9 U', '+12 U'))
        plt.yticks(fontsize=36, fontweight='bold')

        ax1 = inset_axes(ax, width="15%",  height="3%", loc='lower right')
        cbar = plt.colorbar(scatter_tal, cax=ax1, orientation="horizontal")
        ax1.xaxis.set_ticks_position("top")
        ax1.xaxis.set_tick_params(width=5)

    if graph_type == 'bars':

        covi = [0.0 for x in range(len(es_read_matrix[0][l_from:l_to]))]
        edii = [0.0 for x in range(len(covi))]
        xi   = [10.0*x for x in range(len(covi))]
        for y in range(len(es_read_matrix)):
            for x in range(len(covi)):
                if y == 20:
                    covi[x] += es_read_matrix[y][x]
                else:
                    edii[x] += es_read_matrix[y][x]

        font = {'family': 'serif', 'color':  'black', 'weight': 'normal', 'size': 56 }
        p1 = plt.bar(xi, covi, width=10.0, alpha=1)
        p2 = plt.bar(xi, edii,bottom=covi, width=10.0, alpha=1)
        ax = plt.gca()
        x_ticks_list = [x for x in range(10 * 25, 10 * ref_tless_length - 100, 10 * 25)]
        x_ticks_labels = ['{}'.format(int(x/10.0)) for x in x_ticks_list]; x_ticks_labels[-1] = x_ticks_labels[-1] + ' nonT'
        ax.set_xticks(x_ticks_list)
        ax.set_xticklabels(x_ticks_labels, fontdict=font)
        plt.xlim(0, 10.0 * ref_tless_length)
        ax.tick_params(axis='y', which='both', labelsize=56)
        ax.set_ylabel('Reads',fontdict=font)

        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        legend_elements = [matplotlib.lines.Line2D([0], [0], color=colors[0], lw=60, label='unedited reads'),
                   matplotlib.lines.Line2D([0], [0], color=colors[1], lw=60, label='edited reads')]

        ax.legend(handles=legend_elements, loc='upper right', fontsize=56)


    output_pic = 'output.pdf'
    if runArgs.o is not None:
        output_pic = str(runArgs.o)
    plt.savefig(output_pic)
