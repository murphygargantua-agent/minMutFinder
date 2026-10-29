#!/usr/bin/env python3

import gzip
import sys
import os
import shutil
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from read_count import read_count


def _open_fastq(path):
    """Open a FASTQ file, handling both gzipped and plain text."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


for i in range(len(sys.argv)):
    if sys.argv[i] == '--out-dir':
        out_dir = sys.argv[i + 1]
    elif sys.argv[i] == '--r1':
        r1 = sys.argv[i + 1]
    elif sys.argv[i] == '--r2':
        r2 = sys.argv[i + 1]
    elif sys.argv[i] == '--ref-seq':
        ref_seq = sys.argv[i + 1]
    elif sys.argv[i] == '--sample':
        SAMPLE = sys.argv[i + 1]
    elif sys.argv[i] == '--r1qc':
        FQ_R1 = sys.argv[i + 1]
    elif sys.argv[i] == '--r2qc':
        FQ_R2 = sys.argv[i + 1]
    elif sys.argv[i] == '--r1UPqc':
        FQ_R1_UP = sys.argv[i + 1]
    elif sys.argv[i] == '--r2UPqc':
        FQ_R2_UP = sys.argv[i + 1]


FASTQ = out_dir + '/fastq'
QC_DIR = out_dir + '/qc'
qc_metrics_file = QC_DIR + '/qc_metrics.csv'
PLOTS = out_dir + '/plots'
ASSEMBLY = out_dir + '/assembly'


R1_file, R2_file = SeqIO.parse(_open_fastq(r1), 'fastq'), SeqIO.parse(_open_fastq(r2), 'fastq')
R1_qc_paired_file = SeqIO.parse(_open_fastq(FQ_R1), 'fastq')
R2_qc_paired_file = SeqIO.parse(_open_fastq(FQ_R2), 'fastq')
R1_qc_file = SeqIO.parse(_open_fastq(FQ_R1_UP), 'fastq')
R2_qc_file = SeqIO.parse(_open_fastq(FQ_R2_UP), 'fastq')

reads_total = read_count(r1, r2)
reads_filt_qc_unpaired = read_count(FQ_R1_UP, FQ_R2_UP)
reads_filt_qc_paired = read_count(FQ_R1, FQ_R2)
reads_filt_qc = int(reads_filt_qc_paired) + int(reads_filt_qc_unpaired)

qc_metrics = None
if not os.path.isfile(QC_DIR + '/qc_metrics.csv'):
    qc_metrics = open(QC_DIR + '/qc_metrics.csv', 'w')
    print('QC metrics file created')
    qc_metrics.write('sample,test,score' + '\n')
    qc_metrics.write(SAMPLE + ';reads_total;' + str(reads_total) + '\n')
    qc_metrics.write(SAMPLE + ';reads_filt_qc;' + str(reads_filt_qc) + '\n')
    qc_metrics.write(SAMPLE + ';reads_filt_qc_paired;' + str(reads_filt_qc_paired) + '\n')
if qc_metrics is not None:
    qc_metrics.close()

qc = pd.read_csv(qc_metrics_file, sep=',')

sns.set_style("darkgrid", {"axes.facecolor": ".9"})
sns.set_context("paper")
pl = sns.barplot(x=qc.test, y=qc.score)
pl.tick_params(axis='x', labelsize=5)
plt.savefig(PLOTS + '/qc_metrics.svg', format='svg', dpi=2400)