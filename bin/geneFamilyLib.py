from pycbio.hgdata.bed import Bed
from collections import defaultdict, deque

def whackVersion(accid):
    return accid.split('.')[0]

class GeneBoundsBed(Bed):
    """adds properties that are saved in extraCols record
    geneSym - symbol (HGNC)
    hgncId - symbol id
    geneIds - database gene ids, there maybe multiple mapping to the same symbol with readthroughs.
    geneType - type of gene
    """

    def __init__(self, chrom, chromStart, chromEnd, name=None, score=None, strand=None,
                 thickStart=None, thickEnd=None, itemRgb=None, blocks=None, extraCols=None):
        super().__init__(chrom, chromStart, chromEnd, name=name, score=score, strand=strand,
                         thickStart=thickStart, thickEnd=thickEnd, itemRgb=itemRgb, blocks=blocks, extraCols=extraCols)
        self._geneIdsCache = None
        if self.extraCols is None:
            self.extraCols = tuple(4 * [""])

    @property
    def geneSym(self):
        return self.extraCols[0]

    @property
    def hgncId(self, v):
        self.extraCols[1] = v

    def _mkGeneIdsCache(self):
        self._geneIdsCache = tuple(self.extraCols[2].split(','))

    @property
    def geneIds(self):
        if self._geneIdsCache is None:
            self._mkGeneIdsCache()
        return self._geneIdsCache

    def addGeneIds(self, newIds):
        if self._geneIdsCache is None:
            self._mkGeneIdsCache()
        self._geneIdsCache = tuple(set(self._geneIdsCache + tuple(newIds)))
        ex = list(self.extraCols)
        ex[2] = ",".join(sorted(self._geneIdsCache))
        self.extraCols = tuple(ex)

    def addGeneId(self, newId):
        self.addGeneIds((newId, ))

    @property
    def geneType(self):
        return self.extraCols[3]

    @classmethod
    def parse(cls, row, numStdCols=None):
        assert numStdCols == 6
        return Bed.parse(row, numStdCols=6)


# code to merge into gene bounds
def overlaps(b1, b2):
    return (b1.chrom == b2.chrom) and (b1.strand == b2.strand) and (b1.start < b2.end) and (b1.end > b2.start)

def merge(b1, b2):
    b = GeneBoundsBed(b1.chrom, min(b1.start, b2.start), max(b1.end, b2.end), b1.name, 0, b1.strand, extraCols=b1.extraCols)
    b.addGeneIds(b2.geneIds)
    return b

def mergeOrAdd(outRanges, inBed):
    for i in range(len(outRanges)):
        if overlaps(outRanges[i], inBed):
            outRanges[i] = merge(outRanges[i], inBed)
            return True
    outRanges.append(inBed)
    return False

def clusterPass(inRanges, outRanges):
    mergedSome = False
    while len(inRanges) > 0:
        if mergeOrAdd(outRanges, inRanges.pop()):
            mergedSome = True
    return mergedSome

def clusterGeneTrans(transBeds):
    "cluster transcripts for a gene into overlapping beds"
    inRanges = deque(transBeds)
    outRanges = deque()
    while True:
        if not clusterPass(inRanges, outRanges):
            break
        inRanges, outRanges = outRanges, inRanges
    assert len(inRanges) == 0
    return list(outRanges)

def clusterGenes(geneTransBeds):
    "cluster all genes ranges on a chrom so we can sort"
    genesRanges = []
    for transBeds in geneTransBeds.values():
        genesRanges.extend(clusterGeneTrans(transBeds))
    return genesRanges

def setLociNumbers(genesRanges):
    geneIdLociCounts = defaultdict(int)
    for gr in genesRanges:
        gr.score = geneIdLociCounts[gr.name]
        geneIdLociCounts[gr.name] += 1

def buildGeneBounds(geneTransBeds, bedFh):
    genesRanges = clusterGenes(geneTransBeds)
    genesRanges.sort(key=lambda b: (b.chrom, b.start, -b.end))
    setLociNumbers(genesRanges)
    for geneRange in genesRanges:
        geneRange.write(bedFh)
