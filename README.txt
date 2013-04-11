Data analysis pipeline for L1-seq protocol (Ewing and Kazazian 2010)

Prerequisites:
python 2.7+
samtools (http://samtools.sourceforge.net/)
pysam (http://code.google.com/p/pysam/)
tabix (http://samtools.sourceforge.net/tabix.shtml)
bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

Additionally, a bowtie2 index and a fasta index (via samtools faidx) should be
built from the same initial reference genome (in fasta format)

Order of operations:
L1-seq --> FASTQ --> BAM --> per-sample results --> cross-sample results --> primer design

Starting with: 
1. FASTQs (example.fastq.gz)
2. a bowtie2 reference (/path/to/bowtie2ref/hg19)
3. an indexed fasta reference (/path/to/fasta/hg19.fa.fai)

1. make BAM file (output will be example.bam):
./run_bowtie2.py example.fastq.gz /path/to/bowtie2ref/hg19 /path/to/fasta/hg19.fa.fai

2. analyze BAM file (l1seq.py):
usage: l1seq.py [-h] -b INBAMFILE [--minaln MINALEN] [--minmapq MINMAPQ]
                [--readwindow READWINDOW]

optional arguments:
  -h, --help                     show this help message and exit
  -b INBAMFILE, --bam INBAMFILE  input BAM from L1-seq experiment
  --minaln MINALEN               minimum alignment length (default=50)
  --minmapq MINMAPQ              minimum mapq (default=10)
  --readwindow READWINDOW        maximum distance between reads in clusters (default=300bp)


To use defaults (recommended) and save output to example.l1seq.txt:
./l1seq.py example.bam > example.l1seq.txt

3. index the results:
bgzip example_result.l1seq.txt
tabix -s 1 -b 3 -e 4 example_result.l1seq.txt.gz

4. compare results (and save in example_comparison.tsv)
now say there are 3 results: 
example1.l1seq.txt.gz
example2.l1seq.txt.gz
example3.l1seq.txt.gz

./compare.py example1.l1seq.txt.gz example2.l1seq.txt.gz example3.l1seq.txt.gz > example_comparison.tsv


