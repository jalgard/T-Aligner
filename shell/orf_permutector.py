'''
    Author: jalgard
'''

import sys
import argparse
import matplotlib
import os
import itertools
from collections import defaultdict
from suffix_tree import Tree


def log(what, out_dest = sys.stderr.writelines, logger_name = "Permutector", stage = ""):
    """ Fancy logging function """

    out_dest('{} --[{}]--> {}\n'.format(logger_name, stage, what))


def get_opts(args):
    """ Parse input arguments """

    parser = argparse.ArgumentParser()

    parser.add_argument('--r', help="Cryptogene reference sequence files")
    parser.add_argument('--t', help="Read alignments produced by 'alignlib' tool, TAF format")
    parser.add_argument('--o', help='Predicted transcripts in FASTA format')

    return parser.parse_args(args)


def read_ref(file):
    """ Read single-entry FASTA file extracting sequence and name """

    single_fa_file = open(file, 'r')

    data = single_fa_file.readlines()
    if len(data) < 2 or data[0][0] != '>' or '>' in ''.join(data[1:]):
        log('Wrong or empty reference file', stage = 'FASTA READER')
        exit()

    single_fa_file.close()

    return ''.join(line.strip() for line in data[1:]).upper(), data[0].strip()[1:]


def get_crypto_tl(crypto_seq):
    """ Count cryptogene T-less length """

    return len(crypto_seq) - crypto_seq.count('T')


def read_taf_to_matrix(file, crypto_tl):
    """ Read TAF read alignments into reference matrix """

    log('Processing mapped reads edits', stage='TAF-READER')
    taf_matrix = [defaultdict(lambda : 0) for x in range(crypto_tl)]

    with open(file, 'r') as taf:
        for line in taf:
            taf_tokens = line.strip().split('\t')
            aln_start = int(taf_tokens[2])

            for p, ess in enumerate(taf_tokens[4][:-1].split(';')):
                taf_matrix[aln_start + p]['{}{}'.format('+' if int(ess) >0 else '-', abs(int(ess))*'U')] += 1

    return taf_matrix


def read_taf_to_suffix_tree(file):
    """ Preprocess mapped reads to suffix tree """

    log('Building suffix tree from reads data', stage='TAF-READER')
    suffix_tree = Tree()
    rc = 0

    with open(file, 'r') as taf:
        for line in taf:
            rc += 1
            taf_tokens = line.strip().split('\t')
            suffix_tree.add(str(rc), taf_tokens[5])

    return suffix_tree


def make_tl_dna(crypto_seq):
    """ Generate T-less reference profile for cryptogene """

    t_acc = 0
    tlp = []
    tls = crypto_seq.replace('T', '')

    while crypto_seq[0] == 'T':
        crypto_seq = crypto_seq[1:]

    for nuc in crypto_seq[1:]:
        if nuc == 'T':
            t_acc += 1

        else:
            tlp.append(t_acc)
            t_acc = 0

    tlp.append(t_acc)

    return tlp, tls


def only_print_naive_permut_input(crypto_seq, taf_editing_matrix):
    """ Prints list of lists for naive permutations with support values """

    permut_input = []
    tlp, tls = make_tl_dna(crypto_seq)

    for nuc, ref_edit, edits in zip(tls, tlp, taf_editing_matrix):
        states = []

        for edit in edits:
            if edit == '-':
                states.append('{}{} |{}'.format(nuc, 'T' * ref_edit, edits[edit]))

            elif edit[0] == '-' and edits[edit] > 4:
                 states.append('{}{} |{}'.format(nuc, 'T' * ( ref_edit - len(edit) + 1 ), edits[edit]))

            elif edit[0] == '+' and edits[edit] > 4:
                 states.append('{}{} |{}'.format(nuc, 'T' * ( ref_edit + len(edit) - 1 ), edits[edit]))
        states[1:] = sorted(states[1:], key=lambda x : int(x.split('|')[1]), reverse=True)

        permut_input.append(states)

    print(permut_input)


def make_naive_permut_input(crypto_seq, taf_editing_matrix, supp_threshold = 3):
    """ Generate a list of lists for permutation with naive algorithm """

    permut_input = []
    tlp, tls = make_tl_dna(crypto_seq)

    for nuc, ref_edit, edits in zip(tls, tlp, taf_editing_matrix):
        states = []

        for edit in edits:
            if edit == '-':
                states.append('{}{}'.format(nuc, 'T' * ref_edit))

            elif edit[0] == '-' and edits[edit] > supp_threshold:
                 states.append('{}{}'.format(nuc, 'T' * ( ref_edit - len(edit) + 1 )))

            elif edit[0] == '+' and edits[edit] > supp_threshold:
                 states.append('{}{}'.format(nuc, 'T' * ( ref_edit + len(edit) - 1 )))
        permut_input.append(states)

    return permut_input


def compress_permut_input(permut_input):
    """ Joins lists into longer strings if no alternative path exists """

    compressed = [permut_input[0]]
    i = 1
    while i < len(permut_input):
        if len(permut_input[i]) == 1:
            compressed.append(permut_input[i])
            i += 1
            while i < len(permut_input) and len(permut_input[i]) == 1:
                compressed[-1][0] += permut_input[i][0]
                i += 1
            compressed.append(permut_input[i])

        else:
            compressed.append(permut_input[i])
        i += 1

    return compressed


def filter_permutation_list(permutations, support, checked = 0):
    """ Returns the list of only supported permutations of the sequence """

    filtered = []
    log('Filtering permutations')
    for p in permutations:
        if len(p) < 50:
            if support.find(p):
                filtered.append(p)

        else:
            supported = True
            for x in range(checked, len(p)-40, 40):
                if support.find(p[x:x + 40]) == False:
                    supported = False
                    break

            if supported:
                filtered.append(p)

    log('Filtering permutations: {} excluded from total {}'.format(len(permutations) - len(filtered), len(permutations)))
    return filtered

'''
def stepwise_permutations(permut_input, step_size = 3):
    """ Compute permutations in bins of fixed size """

    stepwise_results = []

    for i in range(0, len(permut_input), step_size):
        step_result = list(itertools.product(*permut_input[i:i+step_size]))
        seqs = [''.join(seq) for seq in step_result]
        stepwise_results.append(seqs)

    return stepwise_results


def hierarchical_stepwise_permutations_loop(permut_input, support, step_size = 8):
    """ Do bottom-top joining of stepwise permutations """

    results = permut_input
    log('Starting hierarchical permutation process')
    while len(results) >= step_size:
        log('Next phase with size {}'.format(len(results)))
        results = stepwise_permutations(results, step_size)
        filtered_results = []
        for r in results:
            filtered_results.append(filter_permutation_list(r, support))

        results = filtered_results
        step_size = 2

    return results
'''

def step_by_step_walk(compressed_permut_input, support):
    """ Do step-by-step movement through the input checking growing ORFs each step """

    step = 1
    checked = 0
    current_seeds = compressed_permut_input[0]
    while step < len(compressed_permut_input):
        step_result = list(itertools.product(*[current_seeds,compressed_permut_input[step]]))
        current_seeds = filter_permutation_list([''.join(seq) for seq in step_result], support, checked)
        checked = int(len(current_seeds[0]) / 1.7)
        step += 1

    return current_seeds


def permute_orfs():
    """ Run main alogrithm of ORFs generation via permutations """

    args = get_opts(sys.argv[1:])

    if args.r is None or os.path.isfile(args.r) is False:
        log('Error reading reference file or --r option not provided', stage = 'MAIN')

    crypto_seq, ref_name = read_ref(args.r)
    crypto_tl = get_crypto_tl(crypto_seq)

    if args.t is None or os.path.isfile(args.t) is False:
        log('Error reading alignments file or --t option not provided', stage = 'MAIN')

    taf_editing_matrix = read_taf_to_matrix(args.t, crypto_tl)

    log('Computing initial permutation input list')
    permut_input = make_naive_permut_input(crypto_seq, taf_editing_matrix)
    log('Compressing permutation input list')
    compressed_permut_input = compress_permut_input(permut_input)
    only_print_naive_permut_input(crypto_seq, taf_editing_matrix)
    #print(compressed_permut_input)
    #exit()

    log('Making suffix tree')
    suff_tree = read_taf_to_suffix_tree(args.t)
    result = step_by_step_walk(compressed_permut_input[:-1], suff_tree)
    with open('results.fasta', 'w') as ofile:
        n = 0
        for r in result:
            n += 1
            ofile.writelines('>{}\n{}\n'.format(n, r))


if __name__ == '__main__':
    """ Main programm flow """

    log('Starting protocol')
    permute_orfs()
