import sys, argparse, matplotlib, os
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np
from matplotlib import gridspec

# Script implements Dr. Sara's request

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

    optionParser.add_argument('--m', '--main',
        action='store',
        help="File with main path (canonical ORF), TAF format")

    optionParser.add_argument('--t', '--taf',
        action='store',
        help="File with read alignments in TAF format")

    optionParser.add_argument('--n', '--ref',
        action='store',
        help="Name of the fasta record")

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

    for y in range(-20, 20):
        es_read_matrix.append([0 for x in range(ref_tless_length)])

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

            for p, ess in enumerate(es):
                if s + p >= ref_tless_length:
                    continue
                es_read_matrix[20-ess][s + p] += 1

    mainMrnaFile = None
    if runArgs.m is not None:
        if os.path.isfile(runArgs.m):
            mainMrnaFile = runArgs.m
    else:
        sys.stderr.write('No TAF for main path!')
        exit()

    fig_font_size = 44
    fig=plt.figure(figsize=[ref_tless_length / 7, 11])
    plt.gcf().set_dpi(100)
    spec = gridspec.GridSpec(ncols=1, nrows=3, height_ratios=[1, 2, 5])
    ax_main = fig.add_subplot(spec[2])
    ax_prod = fig.add_subplot(spec[1], sharex = ax_main)
    ax_prec = fig.add_subplot(spec[0], sharex = ax_main)
    ax_prec.set_ylim(0, 100)

    ax_prec.xaxis.set_major_locator(plt.NullLocator())
    ax_prec.xaxis.set_major_formatter(plt.NullFormatter())

    covi = [0.0 for x in range(len(es_read_matrix[0]))]
    edii = [0.0 for x in range(len(covi))]
    xi   = [10.0*x for x in range(len(covi))]

    main_orf_es = ['-' for x in range(len(covi))]
    with open(mainMrnaFile, 'r') as ifile:
        for line in ifile:
            toks = line.split('\t')
            s = int(toks[2])
            es = [int(x) for x in toks[4][:-1].split(';')[1:-2]]
            for p, ess in enumerate(es):
                main_orf_es[s + p] = 20-ess

    not_edited_main = [0.0 for x in range(len(covi))]

    for y in range(len(es_read_matrix)):
        for x in range(len(covi)):
            if main_orf_es[x] == '-' or main_orf_es[x] == 20:
                covi[x] = 0
                edii[x] = 0
                not_edited_main[x] = 100.0
            elif y == main_orf_es[x]:
                edii[x] += es_read_matrix[y][x]
            elif y != 20:
                covi[x] += es_read_matrix[y][x]


    precent_prod = [float(100.0 * edii[i] / (1.0+edii[i] + covi[i])) for i in range(len(covi))]
    max_prod = max([x+y for x,y in zip(covi, edii)])
    ax_prod.set_ylim(0, max_prod)
    not_edited_main_prod = [0.0 if i == 0.0 else max_prod for i in not_edited_main]

    p1 = ax_prod.bar(xi, covi, color='#8e0152', width=10.0, alpha=1)
    p2 = ax_prod.bar(xi, edii, color='#de77ae', bottom=covi, width=10.0, alpha=1)
    p3 = ax_prod.bar(xi, not_edited_main_prod, color='#bababa', width=10.0, alpha=1)
    x_ticks_list = [x for x in range(10 * 25, 10 * ref_tless_length - 100, 10 * 25)]
    x_ticks_labels = ['{}'.format(int(x/10.0)) for x in x_ticks_list]; x_ticks_labels[-1] = x_ticks_labels[-1] + ' nonT'
    ax_prod.set_xticks(x_ticks_list)
    ax_prod.set_xticklabels(x_ticks_labels)
    plt.xlim(0, 10.0 * ref_tless_length)
    ax_prod.tick_params(axis='y', labelsize=fig_font_size)
    ax_prod.set_ylabel('Editing',fontsize=fig_font_size) #fontweight='bold'
    ax_prod.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax_prec.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

    #legend_elements_prod = [matplotlib.lines.Line2D([0], [0], color='#de77ae', lw=5, label='Canonical'),
    #        matplotlib.lines.Line2D([0], [0], color='#8e0152', lw=5, label='Non-canonical'),
    #        matplotlib.lines.Line2D([0], [0], color='#bababa', lw=8, label='Unedited sites')]
    #ax_prod.legend(loc='upper left', handles=legend_elements_prod, fontsize=fig_font_size)

    covi = [0.0 for x in range(len(es_read_matrix[0]))]
    edii = [0.0 for x in range(len(covi))]
    xi   = [10.0*x for x in range(len(covi))]
    for y in range(len(es_read_matrix)):
        for x in range(len(covi)):
            if y == 20:
                covi[x] += es_read_matrix[y][x]
            else:
                edii[x] += es_read_matrix[y][x]


    p4 = ax_main.bar(xi, covi, color='#2166ac', width=10.0, alpha=1)
    p5 = ax_main.bar(xi, edii, color='#92c5de', bottom=covi, width=10.0, alpha=1)
    ax_main.tick_params(axis='y', labelsize=fig_font_size)
    ax_main.tick_params(axis='x', which='both', labelsize=fig_font_size)
    ax_main.set_ylabel('Coverage',fontsize=fig_font_size)

    legend_elements_main = [matplotlib.lines.Line2D([0], [0], color='#92c5de', lw=15, label='Edited'),
            matplotlib.lines.Line2D([0], [0], color='#2166ac', lw=15, label='Unedited'),
            matplotlib.lines.Line2D([0], [0], color='#de77ae', lw=15, label='Canonical'),
            matplotlib.lines.Line2D([0], [0], color='#8e0152', lw=15, label='Non-canonical'),
            matplotlib.lines.Line2D([0], [0], color='#bababa', lw=15, label='Unedited sites')]
    ax_main.legend(loc='upper left', handles=legend_elements_main, fontsize=fig_font_size, frameon=False)
    p7 = ax_prec.bar(xi, precent_prod, color='#de77ae', width=10.0, alpha=1)
    p8 = ax_prec.bar(xi, not_edited_main, color='#bababa', width=10.0, alpha=1)
    ax_prec.set_ylabel('   %',fontsize=fig_font_size, fontweight='bold')
    ax_prec.tick_params(axis='y', labelsize=fig_font_size)
    #legend_elements_prec = [matplotlib.lines.Line2D([0], [0], color='#de77ae', lw=5, label='% canonical')]
    #ax_prec.legend(loc='upper left', handles=legend_elements_prec, fontsize=fig_font_size)

    fig.tight_layout()

    output_pic = 'output.pdf'
    if runArgs.o is not None:
        output_pic = str(runArgs.o)
    plt.savefig(output_pic)
