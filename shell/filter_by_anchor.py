import sys, os

def trim_anchor(grna, mrna):
	alnstr  = ''
	for g, m in zip(grna, mrna):
		if g != m:
			if g == 'A' and m == 'G':
				alnstr += ':'
			elif g == 'C' and m == 'T':
				alnstr += ':'
			else:
				alnstr += '^'
		else:
			alnstr += '|'

	p = len(alnstr) - 1
	anclen = 0
	t = 0
	while p > 0 and anclen < 4:
		if alnstr[p] == '|':
			anclen += 1
		else:
			anclen = 0
		p -= 1
		t += 1
	return t, alnstr[0:p+anclen+1]

def load_fasta(ff):
	data = {}
	last = ''
	with open(ff, 'r') as ifile:
		for line in ifile:	
			if line[0] == '>':
				last = line[1:].rstrip()
				data[last] = ''
			else:
				data[last] += line.rstrip().upper()
	return data

from collections import defaultdict as dd
gdata = dd(lambda : dd(lambda : []))

miniseq = load_fasta(sys.argv[4])
input_data = sys.argv[1]
min_len_cutoff = int(sys.argv[2])
gene_name = sys.argv[3]


for line in open(input_data, 'r').readlines():
	toks = line.split('\t')
	trim = trim_anchor(toks[5], toks[4])
	if len(trim[1]) < min_len_cutoff:
		continue
	mm = trim[1].count('^')
	gu = trim[1].count(':')
	pstart = miniseq[toks[2]].find(toks[5])
	mrna = toks[4][0:len(trim[1])]
	grna = toks[5][0:len(trim[1])]
	print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(toks[0], toks[1], toks[2], toks[3], mrna, grna, gu, mm, pstart))



