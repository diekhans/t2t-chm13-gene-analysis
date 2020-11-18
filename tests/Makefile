
.SECONDARY:

test: geneBoundsBlatTests geneBoundsGencodeTests \
	genePredFrameShiftsTests catBoundsAddBlatOverTest missingGenesTests


geneBoundsBlatTests: geneBoundsTest geneBoundsSymTest geneBoundsExpandTest geneBoundsMinExonsTest

geneBoundsTest: mkdirs output/gencode.bigBlatPsl
	../bin/geneBoundsFromBigBlatPsl --nameField=geneId --hgncOnly output/gencode.bigBlatPsl output/$@.bed
	diff expected/$@.bed output/$@.bed

geneBoundsSymTest: mkdirs output/gencode.bigBlatPsl
	../bin/geneBoundsFromBigBlatPsl --nameField=geneSym output/gencode.bigBlatPsl output/$@.bed
	diff expected/$@.bed output/$@.bed

geneBoundsExpandTest: mkdirs output/gencode.bigBlatPsl
	../bin/geneBoundsFromBigBlatPsl --sourceAnnot=input/gencode.source.bed --maxSourceExpansion=1.5 output/gencode.bigBlatPsl output/$@.bed
	diff expected/$@.bed output/$@.bed

geneBoundsMinExonsTest: mkdirs output/gencode.bigBlatPsl
	../bin/geneBoundsFromBigBlatPsl --sourceAnnot=input/gencode.source.bed --minExons=3 output/gencode.bigBlatPsl output/$@.bed
	diff expected/$@.bed output/$@.bed

output/gencode.bigBlatPsl: input/gencode.bigBlatPslIn input/t2tChm13_20200727.sizes
	@mkdir -p $(dir $@)
	bedToBigBed -tab -type=bed12+17 input/gencode.bigBlatPslIn input/t2tChm13_20200727.sizes $@

genePredFrameShiftsTests: mkdirs
	../bin/genePredFrameShifts --catMetaTsv=input/fakeCat.gp_info input/fakeCat.gp output/$@.tsv output/$@.gp
	diff expected/$@.tsv output/$@.tsv
	diff expected/$@.gp output/$@.gp

catSourceDb =  ../../../build/t2tChm13_20200727/CAT/databases/GRCh38.db 
catBoundsAddBlatOverTest: mkdirs
	../bin/catBoundsAddBlatOver ${catSourceDb} input/catBoundsAddBlatOver.gp_info input/catBoundsAddBlatOver.blat.overstats input/catBoundsAddBlatOver.bounds.bed output/$@.bed
	diff expected/$@.bed output/$@.bed


geneBoundsGencodeTests: geneBoundsGencodeTest

gencodeDb = ../../build/t2t-chm13-v1/CAT/missing/data/gencodeV35.db
geneBoundsGencodeTest: ${gencodeDb} mkdirs
	../bin/geneBoundsForSource --limit=250 --geneType=protein_coding --geneType=lncRNA ${gencodeDb} output/$@.bed


missingGenesTests: missingGeneSourceContextTest

buildBigBed = ../../../T2T-CHM13-hub/bin/buildBigBed
gencodeSrcBigBed = output/gencode-src.bigBed

missingGeneSourceContextTest: ${gencodeSrcBigBed}
	../bin/getMissingGeneSourceContext ${gencodeSrcBigBed} input/gencode-src-missing.bed output/$@.json
	diff expected/$@.json output/$@.json


${gencodeSrcBigBed}: input/gencode-src.bed mkdirs
	${buildBigBed} --extraIndex=name --twoBit=/hive/data/genomes/hg38/hg38.2bit --as=../etc/geneBoundsBed.as  bed9+4 hg38 $@ input/gencode-src.bed

clean:
	rm -rf output

mkdirs:
	@mkdir -p output
