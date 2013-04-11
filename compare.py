#!/usr/bin/env python

import sys
import pysam

'''
compare.py: compares multiple outputs of l1seq.py
'''

class Site:
    def __init__(self, cluster, samples, window=200):
        c = cluster.strip().split()
        self.window = window
        # initialize position, can update later
        self.chrom  = c[0]
        self.start  = int(c[2])
        self.end    = int(c[3])
        self.strand = c[4]
        self.samples = samples
        self.seensample = {}
        for sample in samples:
            self.seensample[sample] = False

        # track info from contributing clusters
        self.inspos = []
        self.starts = []
        self.ends   = []
        self.counts = []
        self.uniqs  = []
        self.widths = []
        self.maptrk = []
        self.mapqs  = []
        self.refelt = []
        self.inelt  = []
        self.nonref = []
        self.ingene = []
        self.inexon = []

    def addsample(self, tabix):
        assert self.seensample[tabix.filename] == False
        self.seensample[tabix.filename] = True

        minpos = self.start - self.window
        if minpos < 0:
            minpos = 0
        maxpos = self.end + self.window

        bestcluster = None
        maxscore = 0

        inspos = start = end = count = uniq = width = 0
        mapscore = mapq = 0.0
        refelt = inelt = nonref = ingene = inexon = 'NA' 

        clusters = [] 
        for cluster in tabix.fetch(self.chrom, minpos, maxpos):
            clusters.append(cluster)

        if len(clusters) > 0:
            if len(clusters) > 1:
                sys.stderr.write("warning: overlapping clusters in region: " + self.chrom +
                                 ":" + str(minpos) + "-" + str(maxpos) + "\n")
                for cluster in clusters:
                    c = cluster.strip().split()
                    count = int(c[5])
                    uniq = int(c[6])
                    if count*uniq > maxscore:
                        maxscore = count*uniq
                        bestcluster = cluster
            else:
                bestcluster = clusters[0]

            if len(clusters) > 1:
                sys.stderr.write("\tbest cluster score was: " + str(maxscore) + "\n")
            strand = bestcluster.strip().split("\t")[4]
            chrom  = bestcluster.strip().split("\t")[0]

            assert chrom == self.chrom
            if self.strand == strand:
                (chrom, inspos, start, end, strand, count, uniq, width, mapscore, 
                 mapq, refelt, inelt, nonref, ingene, inexon) = bestcluster.strip().split("\t")

        self.inspos.append(int(inspos))
        self.starts.append(int(start))
        self.ends.append(int(end))
        self.counts.append(int(count))
        self.uniqs.append(int(uniq))
        self.widths.append(int(width))
        self.maptrk.append(float(mapscore))
        self.mapqs.append(float(mapq))
        self.refelt.append(refelt)
        self.inelt.append(inelt)
        self.nonref.append(nonref)
        self.ingene.append(ingene)
        self.inexon.append(inexon)

    def __str__(self):
        samples = []
        for count in self.counts:
            if count > 0:
                samples.append("YES")
            else:
                samples.append("NO")
        samplestr = "\t".join(samples)

        persample_count = "\t".join(map(str, self.counts))
        persample_uniq  = "\t".join(map(str, self.uniqs))
        persample_width = "\t".join(map(str, self.widths))

        avgcount = str(avg(nonzero(self.counts)))
        avgwidth = str(avg(nonzero(self.widths)))
        avguniq  = str(avg(nonzero(self.uniqs)))
        maxcount = str(max(nonzero(self.counts)))
        maxwidth = str(max(nonzero(self.widths)))
        maxuniq  = str(max(nonzero(self.uniqs)))

        minstart = str(min(nonzero(self.starts)))
        maxend   = str(max(nonzero(self.ends)))

        bestloc = 0
        strand = ''
        if self.strand == '+':
            bestloc = str(minstart)
            strand = 'plus'
        if self.strand == '-':
            bestloc = str(maxend)
            strand = 'minus'

        avgmapscore = str(avg(nonzero(self.maptrk)))
        avgmapq = str(avg(nonzero(self.mapqs)))

        refelt = ','.join(map(str, set(nonNA(self.refelt))))
        inelt  = ','.join(map(str, set(nonNA(self.inelt))))
        nonref = ','.join(map(str, set(nonNA(self.nonref))))
        ingene = ','.join(map(str, set(nonNA(self.ingene))))
        inexon = ','.join(map(str, set(nonNA(self.inexon))))

        return "\t".join((self.chrom, bestloc, minstart, maxend, strand,
                          avgcount, avguniq, avgwidth, maxcount, maxuniq, maxwidth,
                          avgmapscore, avgmapq, refelt, inelt, nonref, ingene, 
                          inexon, samplestr, persample_count, persample_uniq,
                          persample_width))

    def verify(self):
        # make sure all samples have been considered
        for sample in self.seensample.keys():
            if self.seensample[sample] == False:
                return False
        return True

    def match(self,chrom,start,end):
        minpos = self.start - self.window
        maxpos = self.end + self.window
        if self.chrom == chrom:
            if overlap((minpos,maxpos), (start,end)) > 0:
                return True
        return False


def nonzero(a):
    return filter(lambda x:x>0.0, a)


def nonNA(a):
    return filter(lambda x:x!='NA', a)


def avg(a):
    return float(sum(a))/float(len(a))


def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def checksites(sites, chrom,start,end):
    for site in sites:
        if site.match(chrom,start,end):
            return True
    return False


def printheader(samples):
    output = "\t".join(('chrom','pos','start','end','strand','avgcount','avguniq','avgwidth',
                        'maxcount','maxuniq','maxwidth','mapscore','mapq','refelt','inelt',
                        'nonref','ingene','inexon'))
    output += "\t" + "\t".join(samples)
    output += "\t" + "\t".join(map(lambda x:x+"_count", samples))
    output += "\t" + "\t".join(map(lambda x:x+"_uniq", samples))
    output += "\t" + "\t".join(map(lambda x:x+"_width", samples))
    
    print output


if len(sys.argv) > 2:
    tabix = []
    sites = []
    for filename in sys.argv[1:]:
        tabix.append(pysam.Tabixfile(filename))

    for t in tabix:
        for cluster in t.fetch():
            c = cluster.strip().split()
            chrom = c[0]
            start = int(c[2])
            end   = int(c[3])
            if not checksites(sites,chrom,start,end):
                site = Site(cluster, sys.argv[1:])
                for tt in tabix:
                    site.addsample(tt)
                sites.append(site)

    printheader(sys.argv[1:])
    for site in sites:
        assert site.verify()
        print site

else:
    print "usage:",sys.argv[0],"<l1seq tabix 1> <l1seq tabix 2> ... <l1seq tabix n>"
