import sys
import os.path as ops
from pycbio.sys.symEnum import SymEnum, auto
sys.path.append(ops.expanduser("~/kent/src/hg/makeDb/outside/gencode/lib"))
from gencode import biotypes

class BioCategory(SymEnum):
    # in order of preference
    protein_coding = auto()
    lncRNA = auto()
    pseudoGene = auto()
    otherRNA = auto()

    @classmethod
    def fromBioType(cls, bioType):
        if isinstance(bioType, str):
            bioType = biotypes.BioType(bioType)
        if bioType in biotypes.bioTypesCoding:
            return cls.protein_coding
        elif bioType in (biotypes.BioType.lncRNA, biotypes.BioType.macro_lncRNA, biotypes.BioType.bidirectional_promoter_lncRNA):
            return cls.lncRNA
        elif bioType in biotypes.bioTypesPseudo:
            return cls.pseudoGene
        else:
            return cls.otherRNA

    @classmethod
    def fromCatType(cls, catBioType):
        if catBioType == "unknown_likely_coding":
            return cls.protein_coding
        else:
            return cls.fromBioType(catBioType)
