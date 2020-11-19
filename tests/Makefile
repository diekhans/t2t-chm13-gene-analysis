hg38TwoBit = /hive/data/genomes/hg38/hg38.2bit
gencodeDb = ../../build/t2t-chm13-v1/CAT/missing/data/gencodeV35.db
buildBigBed = ../../../T2T-CHM13-hub/bin/buildBigBed


define hg38BigBedCheck
${buildBigBed} --extraIndex=name --twoBit=${hg38TwoBit} --as=../etc/geneBoundsBed.as bed9+4 hg38 output/$@.bigBed output/$@.bed
endef


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

geneBoundsGencodeTest: ${gencodeDb} mkdirs
	../bin/geneBoundsForSource --limit=250 --geneType=protein_coding --geneType=lncRNA ${gencodeDb} output/$@.bed
	${hg38BigBedCheck}

missingGenesTests: missingGeneSourceContextTest missingGeneSourceContextTracksTest \
	missingGeneSourceContextOnlyTrackTest

gencodeSrcBigBed = output/gencode-src.bigBed

missingGeneSourceContextTest: ${gencodeSrcBigBed}
	../bin/getMissingGeneSourceContext ${gencodeSrcBigBed} input/gencode-src-missing.bed output/$@.json
	diff expected/$@.json output/$@.json

missingGeneSourceContextTracksTest: mkdirs
	../bin/missingGeneSourceContextTracks input/missingGeneSourceContext.json output/$@.bed
	diff expected/$@.bed output/$@.bed
	${hg38BigBedCheck}

missingGeneSourceContextOnlyTrackTest: mkdirs
	../bin/missingGeneSourceContextTracks --recType=CONTEXT input/missingGeneSourceContext.json --labelsTsv=output/$@.tsv output/$@.bed
	diff expected/$@.bed output/$@.bed
	diff expected/$@.tsv output/$@.tsv
	${hg38BigBedCheck}

${gencodeSrcBigBed}: input/gencode-src.bed mkdirs
	${buildBigBed} --extraIndex=name --twoBit=${hg38TwoBit} --as=../etc/geneBoundsBed.as bed9+4 hg38 $@ input/gencode-src.bed

clean:
	rm -rf output

mkdirs:
	@mkdir -p output
