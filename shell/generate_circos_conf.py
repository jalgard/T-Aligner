import sys, argparse
from collections import defaultdict as ddict

def GetOptParser():

    optionParser = argparse.ArgumentParser(description='Generate circos config files from gRNA finding data')

    optionParser.add_argument('--grna', dest="grna", action='store',
        help=".FASTA of minicircle contigs, fow strand")

    optionParser.add_argument('--fow', dest="fow", action='store',
        help="Results of gRNA finding, against fow strand")

    optionParser.add_argument('--rev', dest="rev", action='store',
        help="Results of gRNA finding, against rev strand")

    optionParser.add_argument('--mrna', dest="mrna", action='store',
        help=".FASTA of edited mRNA input")

    optionParser.add_argument('--out', dest="output", action='store',
        help="Prefix for filenames to store output")

    return optionParser

runArgs = GetOptParser().parse_args(sys.argv[1:])


def read_fasta(ffile):
    fa_data = []
    with open(ffile, 'r') as ifile:
        for line in ifile:
            if line[0] == '>':
                fa_data.append([line[1:].rstrip().split(' ')[0],''])
            else:
                fa_data[-1][1] += line.rstrip().upper()
    return fa_data


mrna_fa = read_fasta(runArgs.mrna)
grna_fa = read_fasta(runArgs.grna)

with open('{}.runcircos.links'.format(runArgs.output), 'w') as ofile:
    with open(runArgs.fow, 'r') as ifile:
        for line in ifile:
            toks = line.rstrip().split('\t')
            if len(toks) > 8:
                ofile.writelines('{}\t{}\t{}\t{}\t{}\t{}\tident={}\n'.format(toks[0], toks[1], toks[2], toks[4], toks[5], toks[6], "100.0"))
    with open(runArgs.rev, 'r') as ifile:
        for line in ifile:
            toks = line.rstrip().split('\t')
            if len(toks) > 8:
                l=0
                for i in grna_fa:
                    if toks[0] == i[0]:
                        l = len(i[1])
                ofile.writelines('{}\t{}\t{}\t{}\t{}\t{}\tident={}\n'.format(toks[0], l-int(toks[1]), l-int(toks[2]), toks[4], toks[5], toks[6], "100.0"))

with open('{}.runcircos.karyo'.format(runArgs.output), 'w') as ofile:
    for fa in mrna_fa:
        ofile.writelines('chr\t-\t{}\t{}\t0\t{}\tgreen\n'.format(fa[0], fa[0], len(fa[1])))
    for fa in grna_fa:
        ofile.writelines('chr\t-\t{}\t{}\t0\t{}\tgrey\n'.format(fa[0], fa[0], len(fa[1])))

with open('{}.runcircos.conf'.format(runArgs.output), 'w') as ofile:
    ofile.writelines('karyotype = {}.runcircos.karyo\n'.format(runArgs.output))
    ofile.writelines('chromosomes=')
    for fa in mrna_fa + grna_fa:
        ofile.writelines('{};'.format(fa[0]))
    ofile.writelines('\nchromosomes_display_default = no\nchromosomes_units=1\n')
    ofile.writelines('<ideogram>\n<spacing>\ndefault = 0.005r\n</spacing>\nshow_bands=yes\nfill_bands=yes\n')
    ofile.writelines('radius=0.90r\nthickness=100p\nfill=yes\nstroke_color=dgrey\nstroke_thickness=1p\nshow_label=yes\nlabel_size=30\nlabel_radius=dims(image,radius)-60p\nlabel_parallel=yes\n')
    ofile.writelines('</ideogram>\nshow_ticks=no\nshow_tick_labels=no\n')
    ofile.writelines('<plots>\n</plots>\n<links>\n<link>\nfile={}.runcircos.links\n'.format(runArgs.output))
    ofile.writelines('radius=0.94r\nbezier_radius=0.02r\nribbon=yes\ncolor=red_a3\n<rules>\n')
    ofile.writelines('</rules>\n</link>\n</links>\n')
    ofile.writelines('<image>\nangle_offset*=0\n<<include etc/image.conf>>\n</image>\n')
    ofile.writelines('<zooms>\n<zoom>\nchr={}\nstart=0u\nend={}u\nscale=1000\n</zoom>\n</zooms>\n'.format(mrna_fa[0][0], len(mrna_fa[0][1])))
    ofile.writelines('<<include etc/colors_fonts_patterns.conf>>\n<<include etc/housekeeping.conf>>\n')
