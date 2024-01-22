import sys


lines = open(sys.argv[1], 'r').readlines()
lines_out = []


for i in range(len(lines)):
    if '@@ Downstream' in lines[i]:
        main_line = lines[i-7]
        toks = main_line.split('\t')

        # append empty entity as list with fields:
        # lines (list), gu (int), mm (int), length (int), mRNA start (int) and mRNA stop (int)
        lines_out.append([[], int(toks[8]), int(toks[9]), int(toks[10]), int(toks[5]), int(toks[6])])

        for b in range(-7,1):
            lines_out[-1][0].append(lines[i+b][:-1])


lines_out.sort(key=lambda x : x[4])
lines_out.sort(key=lambda x : x[2])

for hit in lines_out:
    for i in hit[0]:
        print('{}'.format(i))
