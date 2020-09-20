

test: geneBoundsTest geneBoundsSymTest geneBoundsExpandTest geneBoundsMinExonsTest \
	genePredFrameShiftsTests catBoundsAddBlatOverTest
#catClustersAddGoTest


# test out not deterministic:  geneFamilyTest
# not useful: geneFamilyExpansionTest

geneFamilyTest: mkdirs
	../bin/geneFamilyBuildInfo --testMode input/ensembl_multigenefamilies.subset.tab output/$@.json
	diff expected/$@.json output/$@.json

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

geneFamilyExpansionTest: mkdirs
	../bin/geneFamilyExpansion --minExpansionCnt=-100 expected/geneFamilyTest.json expected/geneBoundsTest.bed output/$@.bed

genePredFrameShiftsTests: mkdirs
	../bin/genePredFrameShifts --catMetaTsv=input/fakeCat.gp_info input/fakeCat.gp output/$@.tsv output/$@.gp
	diff expected/$@.tsv output/$@.tsv
	diff expected/$@.gp output/$@.gp

catSourceDb =  ../../../build/t2tChm13_20200727/CAT/databases/GRCh38.db 
catBoundsAddBlatOverTest: mkdirs
	../bin/catBoundsAddBlatOver ${catSourceDb} input/catBoundsAddBlatOver.gp_info input/catBoundsAddBlatOver.blat.overstats input/catBoundsAddBlatOver.bounds.bed output/$@.bed
	diff expected/$@.bed output/$@.bed

catClustersAddGoTest: mkdirs
	../bin/catClustersAddGo expected/catBoundsAddBlatOverTest.tsv output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

clean:
	rm -rf output

mkdirs:
	@mkdir -p output
