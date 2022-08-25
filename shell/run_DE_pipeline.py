import os, sys, random, subprocess, string
from collections import defaultdict
from datetime import datetime


# Globals
t_aligner_dir = '/home/jalgard/ta330'
logFile       = 'protocol.log'
def Execute(execute_that):
    PIPE = subprocess.PIPE
    p = subprocess.Popen(execute_that, shell=True, stdin=PIPE, stdout=PIPE,stderr=subprocess.STDOUT, close_fds=True)
    p.communicate()

def Log(log_that):
    log_file = open(logFile, 'a')
    log_file.writelines(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    log_file.writelines(' ->\t')
    log_file.writelines(log_that + '\n')
    log_file.close()


# align reads on the reference
# make DEST files
# run EdgeR on DEST files
# plot data

# example:
# run_DE_pipeline.py R /path/to/reference.py M /path/to/mrna.taf C /path/to/fastq1 C /path/to/fastq2 E /path/to/fastq1 E /path/to/fastq2

plot_main_orf  = False
taf_main_orf   = ''
ref_fasta_file = ''
C_samples = []
E_samples = []
run_descriptor = 'run'

for arg_i in range(1, len(sys.argv)):
    arg = sys.argv[arg_i]
    if arg == 'R':
        arg_i += 1
        ref_fasta_file  = sys.argv[arg_i]
    if arg == 'M':
        plot_main_orf = True
        arg_i += 1
        taf_main_orf  = sys.argv[arg_i]
    if arg == 'C':
        arg_i += 1
        C_samples.append(sys.argv[arg_i])
    if arg == 'E':
        arg_i += 1
        E_samples.append(sys.argv[arg_i])
    if arg == 'D':
        arg_i += 1
        run_descriptor = sys.argv[arg_i]

current_runid = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
logFile       = './DE_{}_{}/protocol_{}_{}.log'.format(run_descriptor, current_runid,run_descriptor, current_runid)
# alignments
Execute('mkdir DE_{}_{}'.format(run_descriptor, current_runid))
pwd = 'DE_{}_{}'.format(run_descriptor, current_runid)

Log('Starting protocol.')
sample_id = 0
Log('Mapping control samples')
for sample in C_samples:
    sample_id += 1
    cmd_align = 'cd {}; {}/bin/alignlib --in_ref {} --in_lib {} --out_prefix sample_{}_'.format(pwd, t_aligner_dir, ref_fasta_file, sample, sample_id)
    Execute(cmd_align)
    Log(cmd_align)
    cmd_dest = 'cd {}; python {}/shell/taf_to_DEST.py --r {} --t sample_{}_mapped_reads.taf --o sample_{}.dest'.format(pwd, t_aligner_dir, ref_fasta_file, sample_id, sample_id)
    Execute(cmd_dest)
    Log(cmd_dest)

Log('Mapping experiment samples')
for sample in E_samples:
    sample_id += 1
    cmd_align = 'cd {}; {}/bin/alignlib --in_ref {} --in_lib {} --out_prefix sample_{}_'.format(pwd, t_aligner_dir, ref_fasta_file, sample, sample_id)
    Execute(cmd_align)
    Log(cmd_align)
    cmd_dest = 'cd {}; python {}/shell/taf_to_DEST.py --r {} --t sample_{}_mapped_reads.taf --o sample_{}.dest'.format(pwd, t_aligner_dir, ref_fasta_file, sample_id, sample_id)
    Execute(cmd_dest)
    Log(cmd_dest)

Log('Running EdgeR...')
sstring = ' '.join('C' * len(C_samples)) + ' ' + ' '.join('E' * len(E_samples))

cmd_edger = 'cd {}; Rscript {}/shell/invoke_EdgeR_on_DEST.R {} {}_{}.results.txt'.format(pwd, t_aligner_dir, sstring, run_descriptor, current_runid)
Execute(cmd_edger)
Log(cmd_edger)

Log('Plotting...')

pl_m = ''
if plot_main_orf == True:
    pl_m = ' --m {}'.format(taf_main_orf)

cmd_plot = 'cd {}; python {}/shell/plot_DEST.py --r {} --d {}_{}.results.txt --o plot_{}_{} {}'.format(pwd, t_aligner_dir, ref_fasta_file, run_descriptor, current_runid, run_descriptor, current_runid, pl_m)
Execute(cmd_plot)
Log(cmd_plot)
