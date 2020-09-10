from pycbio.hgdata.bed import Bed

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
