table geneBounds
"GeneBounds with metadata (bed9+4)"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "alignment Id"
    uint score;         "Score (0-1000), faction identity * 1000"
    char[1] strand;     "+ or - indicates whether the query aligns to the + or - strand on the reference"
    uint thickStart;   "Start of where display should be thick (start codon)"
    uint thickEnd;     "End of where display should be thick (stop codon)"
    uint reserved;     "Used as itemRgb as of 2004-11-22"
    string geneSym;  "gene symbol"
    string geneId;  "gene id"
    string hgncId;  "HGNC id"
    string geneType; "gene type"
    )

