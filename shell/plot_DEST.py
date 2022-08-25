import sys, argparse, matplotlib, os
from matplotlib import pyplot as plt
from matplotlib import cm
from collections import defaultdict
import numpy as np
#import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# INPUT : - reference sequence (cryptogene)
#         - alignments in TAF format
#         - mRNA of main ORF in FAT format
#         - DEST file in EdgeR format
def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--r', '--fa',
        action='store',
        help="Reference fasta file")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output table name")

    optionParser.add_argument('--d', '--dest',
        action='store',
        help="File DE analysis results [EdgeR format]")

    optionParser.add_argument('--s', '--sep',
        action='store',
        help="Separator [default: ,]")

    optionParser.add_argument('--n', '--n',
        action='store',
        help="Cryptogene name [optional]")

    optionParser.add_argument('--m', '--mrna',
        action='store',
        help="Plot main mRNA path from TAF")


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

    referenceSequence = fastaEntries[inputFastaName][1]
    # get cryptogene reference length
    ref_tless_length = 0

    for i in referenceSequence.lstrip().rstrip().upper():
        if i != 'T':
            ref_tless_length += 1

    # fill matrices with zero values
    pval_matrix = []
    fc_matrix = []
    for y in range(-20, 20):
        pval_matrix.append([0.0 for x in range(ref_tless_length)])
        fc_matrix.append([0.0 for x in range(ref_tless_length)])
    # read TAF
    separator = '\t'
    if runArgs.s is not None:
        if os.path.isfile(runArgs.s):
            separator = runArgs.s

    inputDEFile = ''
    if runArgs.d is not None:
        if os.path.isfile(runArgs.d):
            inputDEFile = runArgs.d
    else:
        sys.stderr.write('No DE file found or -d option missing!')
        exit()

    cy = 0
    cx = 0
    de_data = open(inputDEFile, 'r').readlines()[1:]
    for line in de_data:
        toks = line.split(separator)
        cx = int(toks[0].split(':')[0])
        cy = int(toks[0].split(':')[1])
        pval_matrix[cy][cx] = 0.0
        fc_matrix[cy][cx] = 0.0

        if float(toks[4]) < 0.001:
            pval_matrix[cy][cx] = float(toks[4])
            fc = 2**(-1.0*float(toks[1]))
            if fc < 1.0:
                fc = -1.0 / fc
            fc_matrix[cy][cx] = -1 * fc

    output_pic = 'output.pdf'
    if runArgs.o is not None:
        output_pic = str(runArgs.o)

    plot_mrna = False
    mrna_x = []
    mrna_y = []
    if runArgs.m is not None:
        plot_mrna = True
        begin = 0
        data = []
        with open(runArgs.m, 'r') as mrna_taf:
            for line in mrna_taf:
                toks = line.split('\t')
                begin = int(toks[2]) - 1
                data  = [int(x) for x in toks[4].split(';')[:-2]]
        for i, d in enumerate(data):
            mrna_y.append(210 + 10.0 * d)
            mrna_x.append(10.0 * (begin + i))

    xi = []
    yi = []
    zi = []
    for y in range(len(fc_matrix)):
        for x in range(len(fc_matrix[0])):
            if fc_matrix[y][x] != 0.0:
                xi.append(10.0 * x)
                yi.append(410 - 10.0 * y)
                zi.append(fc_matrix[y][x])

    plt.figure(figsize=[ref_tless_length / 6, 12])
    #sns.set_style("whitegrid")
    #sns.set_context("paper")
    plt.plot([-10.0, 10.0 * ref_tless_length + 10.0],[210,210], c='#a1d76a', lw=3, alpha=1.0)
    plt.xlim(-10.0, 10.0 * ref_tless_length + 10.0)
    plt.plot(mrna_x, mrna_y, c='#a1d76a', lw=10, alpha=0.5)
    ax = plt.gca()
    scatter = plt.scatter(xi,yi, c=zi, cmap='bwr', s=200.0, lw=1.0, alpha=1.0).set_zorder(20)

    ax.grid(linestyle='-', linewidth=3, which='both')
    #               -6   -3    0   +3   +6
    ax.set_yticks((90, 120, 150, 180, 210, 240, 270, 300, 330))
    ax.set_yticklabels(('-12 U', '-9 U', '-6 U', '-3 U', 'ref', '+3 U', '+6 U', '+9 U', '+12 U'))
    plt.yticks(fontsize=36, fontweight='bold')

    x_ticks_list = [x for x in range(10 * 25, 10 * ref_tless_length - 100, 10 * 25)]
    x_ticks_labels = ['{} nt'.format(int(x/10.0)) for x in x_ticks_list]
    ax.set_xticks(x_ticks_list)
    ax.set_xticklabels(x_ticks_labels)
    plt.xticks(fontsize=36, fontweight='bold')

    ax1 = inset_axes(ax, width="15%",  height="3%", loc='lower right')
    plt.clim(-10,10)
    cbar = plt.colorbar(scatter, cax=ax1, orientation="horizontal", ticks=[-7, -2, 0, 2, 7])
    ax1.xaxis.set_ticks_position("top")
    ax1.xaxis.set_tick_params(width=5)
    plt.xticks(fontsize=36, fontweight='bold')
    ax1.set_xticklabels(['7', '3', '', '3', '7'])
    #ax.text(10.0 * ref_tless_length - 150, -30, 'POSITION', family="serif", fontsize=22)
    #ax.text(-40, 340, 'EDITING', family="serif", fontsize=22, rotation=90)
    plt.savefig('lfc_{}'.format(output_pic), dpi=300)


    # skip pvalue
    '''
    xi = []
    yi = []
    zi = []
    for y in range(len(pval_matrix)):
        for x in range(len(pval_matrix[0])):
            xi.append(10.0 * x)
            yi.append(410 - 10.0 * y)
            zi.append(-1.0 *  np.log10(pval_matrix[y][x]))

    plt.figure(figsize=[ref_tless_length / 5, 7])
    plt.plot([-10.0, 10.0 * ref_tless_length + 10.0],[210,210], c='#a1d76a', lw=10, alpha=0.3)
    plt.scatter(xi,yi, c=zi, cmap='autumn_r', s=55.0, lw=0.1, alpha=0.95)
    plt.plot(mrna_x, mrna_y, c='#FF0000', lw=1)
    plt.clim(0.0, 0.001)
    plt.colorbar()
    plt.savefig('pval_{}'.format(output_pic), dpi=300)
    '''
