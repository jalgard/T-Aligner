# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import sys, argparse, matplotlib, os
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np


# INPUT : - reference sequence (in forward orientation)
#         - alignments in TAF v2. format for forward strand
#         - alignments in TAF v2. format for revcom strand
def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--r', '--fa',
        action='store',
        help="Reference fasta file [needs reference name if multifasta is given]")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output pdf/png file name")

    optionParser.add_argument('--c', '--taf2',
        action='store',
        help="File with read alignments in TAF format [revcom strand]")

    optionParser.add_argument('--t', '--taf',
        action='store',
        help="File with read alignments in TAF format [given, forward strand]")

    optionParser.add_argument('--n', '--ref',
        action='store',
        help="Name of the fasta record")

    optionParser.add_argument('--s', '--strandness',
        action='store',
        help="Determine graph strandness [ss, un]")

    optionParser.add_argument('--l', '--loc',
        action='store',
        help="Locus of reference to plot [1-based reference coordinates]")

    optionParser.add_argument('--y', '--lim_y',
        action='store',
        help="Limits on y axis [from , to]")



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
    ref_dna_letters = referenceSequence.lstrip().rstrip().upper()
    ref_tless_length = 0
    ref_full_length  = len(ref_dna_letters)

    for i in ref_dna_letters:
        if i != 'T':
            ref_tless_length += 1

            # read TAF

    inputTafFile = ''
    if runArgs.t is not None:
        if os.path.isfile(runArgs.t):
            inputTafFile = runArgs.t
    else:
        sys.stderr.write('No TAF alignments file found or -t option missing!')
        exit()


    coverage_profile_FW = defaultdict(lambda : 0.0)
    coverage_profile_RC = defaultdict(lambda : 0.0)

    inputTafFile_RC = ''
    if runArgs.c is not None:
        if os.path.isfile(runArgs.c):
            inputTafFile_RC = runArgs.c
    else:
        sys.stderr.write('No TAF alignments file found or -c option missing!')
        exit()

    # note that only v2.0 TAF files with 6, 7 fields are accepted!
    with open(inputTafFile, 'r') as ifile:
        for line in ifile:
            toks = line.rstrip().split('\t')
            for pos in range(int(toks[6]), int(toks[7])):
                coverage_profile_FW[pos] += 1
    with open(inputTafFile_RC, 'r') as ifile:
        for line in ifile:
            toks = line.rstrip().split('\t')
            for pos in range(int(toks[6]), int(toks[7])):
                coverage_profile_RC[ref_full_length - pos] += 1


    # only BARS format is now supported
    strandness = 'ss'
    if runArgs.s is not None:
        strandness = str(runArgs.s)

    l_from = 0
    l_to   = ref_full_length
    if runArgs.l is not None:
        f, t = str(runArgs.l).split(',')
        l_from = int(f)
        l_to   = int(t)


    plt.figure(figsize=[12 + 2 * (l_to-l_from) / 1000, 7])

    ax = plt.gca()
    if strandness == 'ss':

        fw_cov = [coverage_profile_FW[x] for x in range(l_from, l_to)]
        rc_cov = [-1.0*coverage_profile_RC[x] for x in range(l_from, l_to)]
        xi   = [x for x in range(l_from, l_to)]

        p_fw = plt.stackplot(xi, fw_cov, colors=('blue'))
        p_rw = plt.stackplot(xi, rc_cov, colors=('blue'))
        ax.hlines(0, 0, l_to-l_from+100, color='black',ls = '--')

    else:

        fw_cov = [coverage_profile_FW[x] for x in range(l_from, l_to)]
        rc_cov = [coverage_profile_RC[x] for x in range(l_from, l_to)]
        xi   = [x for x in range(l_from, l_to)]

        p_un = plt.stackplot(xi, fw_cov, rc_cov, colors=('blue', 'red'))
        ax.hlines(0, 0, l_to-l_from+100, color='black',ls = '--')

    x_lim_min, x_lim_max, y_lim_min, y_lim_max = plt.axis()

    if runArgs.y is not None:
        f, t = str(runArgs.y).split(',')
        y_lim_min = int(f)
        y_lim_max = int(t)
    plt.ylim(-1*y_lim_min, y_lim_max)

    output_pic = 'output.pdf'
    if runArgs.o is not None:
        output_pic = str(runArgs.o)

    # wired, but need to savefig() before y tick labels will actually populate the ticklabels array!!!!
    plt.savefig(output_pic)
    y_labels = [str(x.get_text()).replace('−', '') for x in ax.yaxis.get_ticklabels()]
    ax.set_yticklabels(y_labels)
    plt.savefig(output_pic)
