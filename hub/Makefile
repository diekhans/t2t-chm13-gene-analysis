
.SECONDARY:

hubDir = ../../hub

hubFiles = $(wildcard *.html *.txt */*.html */*.txt)

hubUrl = https://hgwdev.gi.ucsc.edu/~markd/t2t/gene-annotation/hub/hub.txt

all: ${hubFiles:%=${hubDir}/%}
	hubCheck ${hubUrl}


${hubDir}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@

