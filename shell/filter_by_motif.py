import sys

#motif = "AAAGCCTGT"

#motif = sys.argv[2]

#motif_mask = "TAGGCTTTATT|TAGGCTTTACT|TAGGCTTTGTT|TAGGCTTTGCT|CAGGCTTTATT|CAGGCTTTACT|CAGGCTTTGTT|CAGGCTTTGCT|TAGGCTTTATC|TAGGCTTTACC|TAGGCTTTGTC|TAGGCTTTGCC|CAGGCTTTATC|CAGGCTTTACC|CAGGCTTTGTC|CAGGCTTTGCC|TAGGCTTTTATT|TAGGCTTTTACT|TAGGCTTTTGTT|TAGGCTTTTGCT|CAGGCTTTTATT|CAGGCTTTTACT|CAGGCTTTTGTT|CAGGCTTTTGCT|TAGGCTTTTATC|TAGGCTTTTACC|TAGGCTTTTGTC|TAGGCTTTTGCC|CAGGCTTTTATC|CAGGCTTTTACC|CAGGCTTTTGTC|CAGGCTTTTGCC|GTAGGCTTTT|ATAGGCTTTT|CTAGGCTTTT"

motif_mask = "ATAGGCTTT|GTAGGCTTT|ACAGGCTTT|GCAGGCTTT"

def revcom(dna):
	dna_rc = ''
	for i in dna:
		if i == 'T':
			dna_rc = 'A' + dna_rc
		elif i == 'A':
			dna_rc = 'T' + dna_rc
		elif i == 'G':
			dna_rc = 'C' + dna_rc
		elif i == 'C':
			dna_rc = 'G' + dna_rc
	return dna_rc

motifs = [revcom(x) for x in motif_mask.split('|')]

lines = open(sys.argv[1], 'r').readlines()

for i in range(len(lines)):
	if '@@ Downstream' in lines[i]:
		has_motif = False
		for motif in motifs:
			if motif in lines[i]:
				has_motif = True
		if has_motif == True:
			for b in range(-7,1):
				print(lines[i+b][:-1]) 
