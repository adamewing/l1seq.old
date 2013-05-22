#!/usr/bin/env python

import argparse
import pysam

'''Analyze L1 insertion candidates from a BAM file produced by run_bowtie2.py'''

def seqfilter(read, args):
    ''' read is a pysam alignedread '''
    seq = read.seq
    len = float(read.alen)

    if read.is_secondary: # must be primary alignment
        return False

    if read.mapq < int(args.minmapq): # mapping quality cutoff
        return False

    if read.alen < int(args.minalen): # aligned bases cutoff
        return False

    #if read.is_duplicate: # PCR duplicate cutoff
    #    return False

    if seq.count('N') > 3: # ambiguous base cutoff
        return False

    # low complexity / homopolymer cutoff
    for base in ('A', 'T', 'G', 'C'):
        if float(seq.count(base))/len > 0.75:
            return False

    return True

def annotate_cluster(cluster):
    cl_start = min(cluster.starts)
    cl_end   = max(cluster.starts) + 90

    with open("annotation/names.txt", 'r') as annlist:
        for tabixfile in annlist:
            tabix = pysam.Tabixfile("annotation/" + tabixfile.strip(), 'r')
            annots = []
            for entry in tabix.fetch(cluster.chrom, cl_start, cl_end):
                c = entry.strip().split()
                annots.append(c[3])
            if len(annots) > 0:
                annstring = ','.join(set(annots))
                cluster.annotations.append(annstring)
            else:
                cluster.annotations.append('None')

class Cluster:
    def __init__(self, bam, window):
        self.bam = bam
        self.window = int(window)

        # these are set while the cluster is built
        self.chrom = None
        self.starts_plus  = []
        self.starts_minus = []
        self.mapqs = []
        self.lastpos = None 

        # these are set after the cluster is complete
        self.strand = None
        self.starts = []
        self.uniqstarts = None
        self.flagged = False
        self.refline = None
        self.mapscore = 0.0

        self.annotations = []

    def is_empty(self):
        if self.chrom is None:
            return True
        return False

    def add_read(self, read):
        if self.chrom is not None:
            assert self.chrom == self.bam.getrname(read.tid)
        else:
            self.chrom = self.bam.getrname(read.tid)
        if read.is_reverse:
            self.starts_minus.append(read.pos)
        else:
            self.starts_plus.append(read.pos)
        self.mapqs.append(read.mapq)
        self.lastpos = read.pos

    def is_near(self, read):
        if self.chrom != self.bam.getrname(read.tid):
            return False
        if self.lastpos > read.pos - self.window:
            return True
        return False

    def pick_strand(self):
        ''' the strand of the aligned reads in the cluster is opposite the strand
            (or orientation) of the inserted L1 '''
        if len(self.starts_plus) > len(self.starts_minus):
            self.starts = self.starts_plus
            self.strand = '-'
        else:
            self.starts = self.starts_minus
            self.strand = '+'

        # if there is significant disgreement about the strand, flag the cluster
        if len(self.starts_plus) > 0 and len(self.starts_minus) > 0:
            if float(len(self.starts_plus))/float(len(self.starts_minus)) > 0.25:
                self.flagged = True

    def build_uniq_starts(self):
        if not self.starts:
            self.pick_strand()
        self.uniqstarts = set(self.starts)

    def best_position(self):
        if self.strand is None:
            self.pick_strand()
        if self.strand == '+':
            return min(self.starts)
        else:
            return max(self.starts)

    def width(self):
        return max(self.starts)+90 - min(self.starts)

    def get_mapscore(self):
        tabix = pysam.Tabixfile('annotation/hsMap50bp.bed.gz', 'r')
        minpos = min(self.starts)
        maxpos = max(self.starts) + 90
        mapscores = []
        if self.chrom in tabix.contigs:
            for mapregion in tabix.fetch(self.chrom, minpos, maxpos):
                mapscores.append(float(mapregion.strip().split()[3]))

        if mapscores:
            self.mapscore = sum(mapscores)/float(len(mapscores))
        else:
            self.mapscore = 0

    def refL1(self):
        tabix = pysam.Tabixfile('annotation/hg19.primateL1.txt.gz', 'r')
        search_min = search_max = 0
        if self.strand is None:
            self.pick_strand()
        if self.strand == '+':
            search_min = min(self.starts) - self.window*2
            search_max = min(self.starts) + self.window
        else:
            search_min = max(self.starts) - self.window
            search_max = max(self.starts) + self.window*2
        if self.chrom in tabix.contigs:
            # pick least divergent element in strand agreement 
            min_diverge = 1000
            best_match = None

            for nearbyL1 in tabix.fetch(self.chrom, search_min, search_max):
                c = nearbyL1.strip().split()
                refL1_start  = int(c[1])
                refL1_end    = int(c[2])
                refline_name = c[3]
                diverge      = int(c[4])
                refL1_strand = c[5]
                assert refL1_strand in ('-', '+')

                # clusters are explained by ref L1 if the 3' end partially overlaps the bounds +/- window
                if refL1_strand == self.strand and diverge < min_diverge:
                    if refL1_strand == '-' and search_min < refL1_start and search_max > refL1_start:
                        min_diverge = diverge
                        best_match = refline_name
                    if refL1_strand == '+' and search_min < refL1_end and search_max > refL1_end:
                        min_diverge = diverge
                        best_match = refline_name

        self.refline = best_match

    def avgmapq(self):
        return float(sum(self.mapqs))/float(len(self.mapqs))

    def polya(self):
        # try to find 3' junction with poly-A tail
        pass

    def __str__(self):
        locstr = "\t".join(( self.chrom, str(self.best_position()), str(min(self.starts)), str(max(self.starts)+90), self.strand ))
        clrstr = "\t".join(( str(len(self.starts)), str(len(self.uniqstarts)), str(self.width()), str(self.mapscore), str(self.avgmapq()) ))
        annstr = "\t".join(self.annotations)
        outstr = "\t".join((locstr, clrstr, str(self.refline), annstr))
        return outstr

def main(args):
    bam = pysam.Samfile(args.inbamfile, 'rb')
    window = int(args.readwindow)

    cluster = Cluster(bam, window)
    for read in bam.fetch():
        if seqfilter(read, args):
            if cluster.is_empty(): # first read in cluster
                cluster.add_read(read)
            else: # add next read if in range, else start new cluster
                if cluster.is_near(read):
                    cluster.add_read(read)
                else:
                    # output the old cluster
                    cluster.pick_strand()
                    cluster.build_uniq_starts()
                    cluster.refL1()
                    cluster.get_mapscore()
                    annotate_cluster(cluster)
                    print cluster

                    # start a new cluster
                    cluster = Cluster(bam, window)
                    cluster.add_read(read)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-b', '--bam', dest="inbamfile", required=True, help="input BAM from L1-seq experiment")
    parser.add_argument('--minaln', dest="minalen", default=50, help="minimum alignment length  (default=50)")
    parser.add_argument('--minmapq', dest="minmapq", default=10, help="minimum mapq (default=10)")
    parser.add_argument('--readwindow', dest='readwindow', default=300, help="maximum distance between reads in clusters (default=300)")
    
    args = parser.parse_args()
    main(args)



