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

    optionParser.add_argument('--m', '--mrna',
        action='store',
        help="Comma-separated list of transcripts to plot")

    optionParser.add_argument('--t', '--taf',
        action='store',
        help="File with read alignments in TAF v2 format")

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
    referenceSequence = fastaEntries[0][1]   # assumes single-sequence fasta reference
    # get cryptogene reference length
    ref_tless_length = 0
    for i in referenceSequence.lstrip().rstrip().upper():
        if i != 'T':
            ref_tless_length += 1

    inputTafData = None
    if runArgs.t is not None:
        if os.path.isfile(runArgs.t):
            inputTafData = open(runArgs.t, 'r').readlines()
    else:
        sys.stderr.write('No TAF file found or -t option missing!')
        exit()

    takeEntries = []
    if runArgs.m is not None:
        takeEntries = [int(x) for x in (runArgs.m).split(',')]
    else:
        sys.stderr.write('Nothing to plot!')
        exit()


    entry = 0
    tafs = []
    taf_ids = []
    for line in inputTafData:
        entry += 1
        if entry in takeEntries:
            tafs.append(line.rstrip().split('\t'))
            taf_ids.append(entry)

    aligned_data = {'ref' : ''}
    for i in taf_ids:
        aligned_data['mrna_{}'.format(i)] = ''


    max_t = []
    for i in range(ref_tless_length):
        max_t.append(0)

    t_accum = 0
    tlpos   = 0
    for i in referenceSequence.lstrip().rstrip().upper():
        if i != 'T':
            if max_t[tlpos] < t_accum:
                max_t[tlpos] = t_accum
            tlpos += 1
            t_accum = 0
        else:
            t_accum += 1

    for taf in tafs:
        tlpos = int(taf[2])
        t_accum = 0
        for i in taf[5]:
            if i != 'T':
                if max_t[tlpos] < t_accum:
                    max_t[tlpos] = t_accum
                tlpos += 1
                t_accum = 0
            else:
                t_accum += 1

    tlpos = 0
    t_accum = 0
    for i in referenceSequence.lstrip().rstrip().upper():
        if i != 'T':
            aligned_data['ref'] += 'T' * t_accum
            if max_t[tlpos] > t_accum:
                aligned_data['ref'] += '-' * (max_t[tlpos]-t_accum)
            tlpos += 1
            t_accum = 0
            aligned_data['ref'] += i
        else:
            t_accum += 1

    for taf, id in zip(tafs, taf_ids):
        tlpos = int(taf[2])
        for i in range(tlpos):
            aligned_data['mrna_{}'.format(id)] += '-'
            aligned_data['mrna_{}'.format(id)] += '-'*max_t[i]
        t_accum = 0
        for i in taf[5]:
            if i != 'T':
                aligned_data['mrna_{}'.format(id)] += 'T' * t_accum
                if max_t[tlpos] > t_accum:
                    aligned_data['mrna_{}'.format(id)] += '-' * (max_t[tlpos]-t_accum)
                tlpos += 1
                t_accum = 0
                aligned_data['mrna_{}'.format(id)] += i
            else:
                t_accum += 1

    output = sys.stdout
    if runArgs.o is not None:
        output = open(runArgs.o, 'w')

    for aln in aligned_data:
        output.writelines('>{}\n{}\n'.format(aln, aligned_data[aln]))
