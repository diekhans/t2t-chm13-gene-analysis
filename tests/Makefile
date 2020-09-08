

test: geneFamilyTest geneBoundsTest geneBoundsSymTest

# not useful: geneFamilyExpansionTest

geneFamilyTest: mkdirs
	../bin/geneFamilyBuildInfo --testMode input/ensembl_multigenefamilies.subset.tab output/$@.json
	diff expected/$@.json output/$@.json

geneBoundsTest: mkdirs output/gencode.bigBlatPsi
	../bin/geneBoundsFromBigBlatPsl output/gencode.bigBlatPsi output/$@.bed
	diff expected/$@.bed output/$@.bed

geneBoundsSymTest: mkdirs output/gencode.bigBlatPsi
	../bin/geneBoundsFromBigBlatPsl --nameField=geneSym output/gencode.bigBlatPsi output/$@.bed
	diff expected/$@.bed output/$@.bed

output/gencode.bigBlatPsi: mkdirs input/gencode.bigBlatPsiIn input/t2tChm13_20200727.sizes
	bedToBigBed -tab -type=bed12+17 input/gencode.bigBlatPsiIn input/t2tChm13_20200727.sizes $@

geneFamilyExpansionTest: mkdirs
	../bin/geneFamilyExpansion --minExpansionCnt=-100 expected/geneFamilyTest.json expected/geneBoundsTest.bed output/$@.bed


clean:
	rm -rf output

mkdirs:
	@mkdir -p output
