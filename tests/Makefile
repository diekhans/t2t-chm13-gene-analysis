hg38TwoBit = /hive/data/genomes/hg38/hg38.2bit
gencodeDb = ../../build/t2t-chm13-v1/CAT/missing/data/gencodeV35.db
buildBigBed = ../../../T2T-CHM13-hub/bin/buildBigBed


define hg38BigBedCheck
${buildBigBed} --extraIndex=name --twoBit=${hg38TwoBit} --as=../etc/geneBoundsBed.as bed9+4 hg38 output/$@.bigBed output/$@.bed
endef


.SECONDARY:

test: geneBoundsBlatTests geneBoundsGencodeTests \
	genePredFrameShiftsTests catBoundsAddBlatOverTest missingGenesTests \
	mappingStatsTests


geneBoundsBlatTests: geneBoundsTest geneBoundsSymTest

geneBoundsTest: mkdirs
	../bin/geneBoundsForGencodeBed --hgncOnly ${gencodeDb} input/gencode-blat.bed output/$@.bed
	diff expected/$@.bed output/$@.bed

geneBoundsSymTest: mkdirs
	../bin/geneBoundsForGencodeBed  --nameField=geneSym ${gencodeDb} input/gencode-blat.bed output/$@.bed
	diff expected/$@.bed output/$@.bed

genePredFrameShiftsTests: genePredFrameShiftsTest
genePredFrameShiftsTest: mkdirs
	../bin/genePredFrameShifts --catMetaTsv=input/fakeCat.gp_info input/fakeCat.gp output/$@.tsv output/$@.gp
	diff expected/$@.tsv output/$@.tsv
	diff expected/$@.gp output/$@.gp

catSourceDb =  ../../../build/t2tChm13_20200727/CAT/databases/GRCh38.db 
catBoundsAddBlatOverTest: mkdirs
	../bin/catBoundsAddBlatOver ${catSourceDb} input/catBoundsAddBlatOver.gp_info input/catBoundsAddBlatOver.blat.overstats input/catBoundsAddBlatOver.bounds.bed output/$@.bed
	diff expected/$@.bed output/$@.bed

geneBoundsGencodeTests: geneBoundsGencodeTest geneBoundsGencodeBedTest

geneBoundsGencodeTest: ${gencodeDb} mkdirs
	../bin/geneBoundsForSource --limit=250 --geneType=protein_coding --geneType=lncRNA ${gencodeDb} output/$@.bed
	${hg38BigBedCheck}

geneBoundsGencodeBedTest: ${gencodeDb} mkdirs
	../bin/geneBoundsForGencodeBed --geneType=protein_coding --geneType=lncRNA ${gencodeDb} input/gencode-blat.bed --filteredBed=output/$@.filt.bed output/$@.bed
	diff expected/$@.bed output/$@.bed
	diff expected/$@.filt.bed output/$@.filt.bed

missingGenesTests: missingGeneSourceContextTest missingGeneSourceContextTracksTest \
	missingGeneSourceContextOnlyTrackTest

gencodeGeneSrcBigBed = output/gencode-genes-src.bigBed

missingGeneSourceContextTest: ${gencodeGeneSrcBigBed}
	cut -f 4 expected/geneBoundsGencodeBedTest.bed | sort -u > output/$@.blat.genes-ids
	fgrep -f output/$@.blat.genes-ids input/gencode-genes-src.bed > output/$@.src.bed
	../bin/getMissingGeneSourceContext ${gencodeGeneSrcBigBed} output/$@.src.bed  output/$@.json
	diff expected/$@.json output/$@.json

missingGeneSourceContextTracksTest: mkdirs
	../bin/missingGeneSourceContextTracks expected/missingGeneSourceContextTest.json output/$@.bed
	diff expected/$@.bed output/$@.bed
	${hg38BigBedCheck}

missingGeneSourceContextOnlyTrackTest: mkdirs
	../bin/missingGeneSourceContextTracks --recType=CONTEXT expected/missingGeneSourceContextTest.json --labelsTsv=output/$@.tsv output/$@.bed
	diff expected/$@.bed output/$@.bed
	diff expected/$@.tsv output/$@.tsv
	${hg38BigBedCheck}

##
mappingStatsTests: bedMappingStatsTrans bedMappingStatsGene

bedMappingStatsTrans: mkdirs
	../bin/bedMappingStats input/gencode-trans-src.bed input/gencode-blat.bed output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

bedMappingStatsGene: mkdirs
	../bin/bedMappingStats --geneSpan input/gencode-genes-src.bed expected/geneBoundsTest.bed output/$@.tsv
	diff expected/$@.tsv output/$@.tsv


##
${gencodeGeneSrcBigBed}: input/gencode-genes-src.bed mkdirs
	${buildBigBed} --extraIndex=name --twoBit=${hg38TwoBit} --as=../etc/geneBoundsBed.as bed9+4 hg38 $@ input/gencode-genes-src.bed

clean:
	rm -rf output

mkdirs:
	@mkdir -p output
