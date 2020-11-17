from collections import defaultdict, deque
from pycbio.sys.symEnum import SymEnum
from pycbio.hgdata.bed import Bed
from pycbio.hgdata.bed import BedReader
from pycbio.sys import typeOps

def whackVersion(accid):
    return accid.split('.')[0]

def splitStrList(v):
    return v.split(',')

def joinStrList(l):
    return ",".join(sorted(set(l)))

def noneIfNullLike(v):
    return None if v in ("", "NULL", "N/A") else v

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
        if self.extraCols is None:
            self.extraCols = 4 * [""]
        else:
            self.extraCols = ["" if v == "N/A" else v for v in self.extraCols]
        if isinstance(self.extraCols, tuple):
            self.extraCols = list(self.extraCols)

    @property
    def geneSym(self):
        return self.extraCols[0]

    def setGeneSym(self, v):
        self.extraCols[0] = v

    @property
    def hgncId(self):
        return self.extraCols[1]

    def setHgncId(self, v):
        self.extraCols[1] = v

    @property
    def geneIds(self):
        return splitStrList(self.extraCols[2])

    def setGeneIds(self, geneIds):
        self.extraCols[2] = joinStrList(geneIds)

    def addGeneIds(self, newIds):
        self.setGeneIds(set(self.geneIds) | set(newIds))

    def addGeneId(self, newId):
        self.addGeneIds((newId, ))

    @property
    def geneType(self):
        return self.extraCols[3]

    def setGeneType(self, v):
        self.extraCols[3] = v

    @classmethod
    def parse(cls, row, numStdCols=None):
        assert numStdCols == 6
        return super().parse(row, numStdCols=6)

def GeneBoundsBedReader(geneBoundsBed):
    return BedReader(geneBoundsBed, bedClass=GeneBoundsBed, numStdCols=6)

def geneBoundsBigBedRead(bigBedFh, chrom, start, end):
    if start >= end:
        return iter(())
    for entry in typeOps.mkiter(bigBedFh.entries(chrom, start, end)):
        data = entry[2].split('\t')
        yield GeneBoundsBed(chrom, entry[0], entry[1], name=data[0],
                            score=data[1], strand=data[2], extraCols=data[3:])

# code to merge into gene bounds based on overlap
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
    """cluster based on name and overlap, overlapped needed to handle duplication
    geneRanges is a dict by name
    """
    genesRanges = clusterGenes(geneTransBeds)
    genesRanges.sort(key=lambda b: (b.chrom, b.start, -b.end))
    setLociNumbers(genesRanges)
    for geneRange in genesRanges:
        geneRange.write(bedFh)

class NameColumn(SymEnum):
    geneId = 1
    geneSym = 2

def geneBoundsAddCmdOpts(parser):
    parser.add_argument('--nameField', choices=(NameColumn.geneId, NameColumn.geneSym),
                        default=NameColumn.geneId, type=NameColumn,
                        help="which column to use for in name and for grouping into genes")
    parser.add_argument('--hgncOnly', action="store_true")
    parser.add_argument('--geneType', action='append',
                        help='type of gene, maybe repeated')

def geneBoundsProcessCmdOpts(opts):
    if opts.geneType is not None:
        opts.geneType = frozenset(opts.geneType)