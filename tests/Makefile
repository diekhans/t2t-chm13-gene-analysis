

test: geneBoundsTest geneBoundsSymTest geneBoundsExpandTest

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
	../bin/geneBoundsFromBigBlatPsl --sourceAnnot=input/gencode.source.bed output/gencode.bigBlatPsl output/$@.bed
	diff expected/$@.bed output/$@.bed

output/gencode.bigBlatPsl: input/gencode.bigBlatPslIn input/t2tChm13_20200727.sizes
	@mkdir -p $(dir $@)
	bedToBigBed -tab -type=bed12+17 input/gencode.bigBlatPslIn input/t2tChm13_20200727.sizes $@

geneFamilyExpansionTest: mkdirs
	../bin/geneFamilyExpansion --minExpansionCnt=-100 expected/geneFamilyTest.json expected/geneBoundsTest.bed output/$@.bed


clean:
	rm -rf output

mkdirs:
	@mkdir -p output
