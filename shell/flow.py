import sys, argparse, matplotlib, os
from matplotlib import pyplot as plt
from collections import defaultdict
import numpy as np
import networkx as nx


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

    optionParser.add_argument('--t', '--taf',
        action='store',
        help="File with read alignments in TAF format")

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


    nodes_mapping = defaultdict(lambda : {})
    node_counter = 1
    for y in range(-20, 20):
        for x in range(ref_tless_length):
            nodes_mapping[y][x] = node_counter
            node_counter += 1

    FlowGraph = nx.Graph()
    Flow_layout = {}
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
            for p, ess in enumerate(es[0:-1]):
                n_from = nodes_mapping[ess][s+p]
                n_to   = nodes_mapping[es[p+1]][s+p+1]
                if FlowGraph.has_edge(n_from, n_to):
                    FlowGraph[n_from][n_to]['weight'] += 1
                    if ess == 0 and es[p+1] == 0:
                        FlowGraph[n_from][n_to]['weight'] = 1

                else:
                    FlowGraph.add_edge(n_from, n_to, weight=1)

                Flow_layout[n_from] = [10 * (s+p), 205 - 5 * ess]
                Flow_layout[n_to] = [10 * (s+p+1), 205 - 5 * es[p+1]]

    nx.draw_networkx_nodes(FlowGraph, Flow_layout, node_size=1.0, node_color='b')
    floe, few = [], []

    loe, ew = zip(*nx.get_edge_attributes(FlowGraph,'weight').items())
    for x, y in zip(loe, ew):
        if y > np.median(ew):
            floe.append(x)
            few.append(y)
    edata = [[x, y] for x, y in zip(floe, few)]
    edata.sort(key = lambda x: x[1])

    List_of_edges = []
    Edge_weights = []
    for ed in edata:
        List_of_edges.append(ed[0])
        Edge_weights.append(ed[1])
    nx.draw_networkx_edges(FlowGraph, Flow_layout, width=1.0, edgelist=List_of_edges, edge_color=Edge_weights, edge_cmap=plt.cm.Oranges)

    output_pic = 'output.pdf'
    if runArgs.o is not None:
        output_pic = str(runArgs.o)

    plt.savefig(output_pic)
