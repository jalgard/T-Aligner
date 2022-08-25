import sys, argparse, matplotlib, os
import numpy as np
import graphviz
import matplotlib.pyplot as plt
from collections import defaultdict

def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output table name")

    optionParser.add_argument('--p', '--pep',
        action='store',
        help="Fasta with assembled peptides [T-Aligner output file]")

    optionParser.add_argument('--k', '--kmer',
        action='store',
        help="Specify kmer size [8]")

    optionParser.add_argument('--m', '--main',
        action='store',
        help="Specify main orf id [not set]")

    optionParser.add_argument('--i', '--ignore1',
        action='store',
        help="Ignore edges supported by single mRNA")

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

    peptideInputFile = None
    if runArgs.p is not None:
        if os.path.isfile(runArgs.p):
            peptideInputFile = open(runArgs.p, 'r').readlines()
    else:
        sys.stderr.write('Please specify input peptide file')
        exit()

    peptidesFasta = ReadFasta(peptideInputFile)


    graph_nodes = {}
    graph_edges = defaultdict(lambda : 0)
    kms = 8
    if runArgs.k is not None:
        kms = int(runArgs.k)

    peptide_graph = graphviz.Digraph(format='pdf')

    for pep in peptidesFasta:
        frame_scores = {}
        for f in range(kms+1):
            score = 0
            for i in range(f, len(pep[1]) - kms, kms):
                kmer = pep[1][i:i+kms]
                if kmer in graph_nodes:
                    score += 1
            frame_scores[f] = score
        best_frame = max(frame_scores, key=frame_scores.get)
        for i in range(best_frame, len(pep[1]) - kms, kms):
            kmer1 = pep[1][i:i+kms]
            kmer2 = pep[1][i+kms:i+2*kms]
            graph_nodes[kmer1] = 1
            graph_nodes[kmer2] = 1
            graph_edges[kmer1 + ':' + kmer2] += 1


    if runArgs.m is not None:
        main_orf_id = int(runArgs.m)
        main_kmers = {}
        for i in range(0, len(peptidesFasta[main_orf_id][1]) - kms):
            main_kmers[peptidesFasta[main_orf_id][1][i:i+kms]] = 1
        for kmer in main_kmers:
            if kmer in graph_nodes:
                peptide_graph.node(kmer, color = 'lightblue', style='filled')


    ignore1 = False
    if runArgs.i is not None:
        ignore1 = True

    for edge in graph_edges:
        kmer1, kmer2 = edge.split(':')
        if graph_edges[edge] < 2:
            if ignore1 == True:
                continue
            else:
                peptide_graph.edge(kmer1, kmer2)
        elif graph_edges[edge] < 5:
            peptide_graph.edge(kmer1, kmer2, penwidth = '3.0', arrowsize='1.0')
        elif graph_edges[edge] < 10:
            peptide_graph.edge(kmer1, kmer2, penwidth = '4.0', arrowsize='1.0')
        elif graph_edges[edge] < 20:
            peptide_graph.edge(kmer1, kmer2, penwidth = '6.0', arrowsize='1.0')
        elif graph_edges[edge] < 30:
            peptide_graph.edge(kmer1, kmer2, penwidth = '8.0', arrowsize='2.0')
        else:
            peptide_graph.edge(kmer1, kmer2, penwidth = '10.0', arrowsize='2.0')


    output_pic = 'output'
    if runArgs.o is not None:
        output_pic = str(runArgs.o)
    peptide_graph.render(filename=output_pic)
