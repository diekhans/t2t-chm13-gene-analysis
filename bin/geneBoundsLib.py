import json
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

def noneIfEmpty(v):
    return None if v == "" else v

def emptyIfNone(v):
    return "" if v is None else str(v)

class GeneBoundsBed(Bed):
    """adds properties that are saved in extraCols record  (bed9+4)
    geneSym - symbol (HGNC)
    hgncId - symbol id
    geneIds - database gene ids, there maybe multiple mapping to the same symbol with readthroughs.
    geneType - type of gene
    """
    __slots__ = ("geneSym", "hgncId", "geneIds", "geneType")

    def __init__(self, numStdCols, chrom, chromStart, chromEnd, name=None, score=None, strand=None,
                 thickStart=None, thickEnd=None, itemRgb=None,
                 geneSym=None, hgncId=None, geneIds=None, geneType=None,
                 blocks=None, extraCols=None):
        assert numStdCols == 9, "numStdCols expected 9, got {}".format(numStdCols)
        super().__init__(numStdCols, chrom, chromStart, chromEnd, name=name, score=score, strand=strand,
                         thickStart=thickStart, thickEnd=thickEnd, itemRgb=itemRgb,
                         blocks=blocks, extraCols=extraCols)
        if self.itemRgb is None:
            self.itemRgb = ""
        self.geneSym = noneIfEmpty(geneSym)
        self.hgncId = noneIfEmpty(hgncId)
        self.geneIds = []
        if isinstance(geneIds, str):
            self.geneIds.append(geneIds)
        elif geneIds not in (None, ""):
            self.geneIds.extend(geneIds)
        self.geneType = noneIfEmpty(geneType)

    @classmethod
    def create(cls, chrom, chromStart, chromEnd, name, *, score=0, strand=None,
               thickStart=None, thickEnd=None, itemRgb="0,0,0",
               geneSym=None, hgncId=None, geneIds=None, geneType=None):
        return cls(9, chrom, chromStart, chromEnd, name, score, strand,
                   thickStart, thickEnd, itemRgb, geneSym, hgncId, geneIds, geneType)

    @property
    def numColumns(self):
        return self.numStdCols + 4

    def toRow(self):
        return super().toRow() + \
            [emptyIfNone(self.geneSym), emptyIfNone(self.hgncId),
             joinStrList(self.geneIds), emptyIfNone(self.geneType)]

    @classmethod
    def parse(cls, row, numStdCols=None):
        assert numStdCols == 9
        b = super().parse(row[0:9], numStdCols=9)
        b.geneSym = noneIfEmpty(row[9])
        b.hgncId = noneIfEmpty(row[10])
        b.geneIds = splitStrList(row[11])
        b.geneType = noneIfEmpty(row[12])
        return b

def GeneBoundsBedReader(geneBoundsBed):
    return BedReader(geneBoundsBed, bedClass=GeneBoundsBed, numStdCols=9)

def geneBoundsBigBedRead(bigBedFh, chrom, start, end):
    if start >= end:
        return iter(())
    for entry in typeOps.mkiter(bigBedFh.entries(chrom, start, end)):
        data = entry[2].split('\t')
        yield GeneBoundsBed.create(chrom, entry[0], entry[1], data[0], score=data[1], strand=data[2],
                                   thickStart=data[3], thickEnd=data[4], itemRgb=data[5],
                                   geneSym=data[6], hgncId=data[7], geneIds=data[8], geneType=data[9])

class GeneBoundsBedJsonEncoder(json.JSONEncoder):
    jsonCols = ("chrom", "chromStart", "chromEnd", "name",
                "strand", "thickStart", "thickEnd", "itemRgb",
                "geneSym", "hgncId", "geneIds", "geneType")

    def default(self, obj):
        if isinstance(obj, GeneBoundsBed):
            return {k: getattr(obj, k) for k in self.jsonCols}
        return json.JSONEncoder.default(self, obj)


# code to merge into gene bounds based on overlap
def overlaps(b1, b2):
    return (b1.chrom == b2.chrom) and (b1.strand == b2.strand) and (b1.start < b2.end) and (b1.end > b2.start)

def merge(b1, b2):
    geneIds = set(typeOps.mkiter(b1.geneIds)) | set(typeOps.mkiter(b2.geneIds))
    b = GeneBoundsBed.create(b1.chrom, min(b1.start, b2.start), max(b1.end, b2.end), b1.name, strand=b1.strand,
                             geneSym=b1.geneSym, hgncId=b1.hgncId,
                             geneIds=list(geneIds), geneType=b1.geneType)
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
    gencodeGeneId = 3  # normally same as geneId

def geneBoundsAddCmdOpts(parser):
    parser.add_argument('--nameField', choices=(NameColumn),
                        default=NameColumn.geneId, type=NameColumn,
                        help="which column to use for in name and for grouping into genes")
    parser.add_argument('--hgncOnly', action="store_true")
    parser.add_argument('--geneType', action='append',
                        help='type of gene, maybe repeated')

def geneBoundsProcessCmdOpts(opts, *, gencodeIdSeparate=False):
    if opts.geneType is not None:
        opts.geneType = frozenset(opts.geneType)
    if (not gencodeIdSeparate) and (opts.nameField == NameColumn.gencodeGeneId):
        opts.nameField = NameColumn.geneId
