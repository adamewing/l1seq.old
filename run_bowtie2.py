#!/usr/bin/env python

import sys
import subprocess
from re import sub

''' wrapper for bowtie2 with parameters tuned for l1-seq'''

def bwt2(fastq, ref, faidx, threads=1):
    # set up filenames
    prefix = fastq.lower()
    for ending in ('.fastq.gz', '.fastq', '.fq', '.fq.gz'):
        if fastq.endswith(ending):
            prefix = sub(ending, '', fastq)
    samfn  = prefix + ".sam"
    bamfn  = prefix + ".bam"
    sortfn = prefix + ".sorted"

    bwt2args = ['bowtie2', '-x' ,ref, '-5', '10', '--local', '--sensitive', '-p', str(threads), '-U', fastq, '-S', samfn]
    bamargs  = ['samtools', 'view', '-bt', faidx, '-o', bamfn, samfn]
    sortargs = ['samtools', 'sort', '-m', '4000000000', bamfn, sortfn]
    idxargs  = ['samtools', 'index', bamfn]


    print "bowtie2 run:", fastq, prefix, threads
    print "mapping:", bwt2args
    subprocess.call(bwt2args)

    print "convert to bam:", bamargs
    subprocess.call(bamargs)

    print "sort bam:", sortargs
    subprocess.call(sortargs)

    print sortfn + ".bam", "-->", bamfn
    os.remove(bamfn)
    os.rename(sortfn + ".bam", bamfn)

    print "index bam:", idxargs
    subprocess.call(idxargs)

if len(sys.argv) == 4:
    bowtie_index = sys.argv[2]
    fasta_index = sys.argv[3]
    assert fasta_index.endswith('.fai')
    bwt2(sys.argv[1], bowtie_index, fasta_index, threads=4)

else:
    print "usage:",sys.argv[0],"<fastq> <bowtie index> <fasta index (.fai from samtools faidx)>"
