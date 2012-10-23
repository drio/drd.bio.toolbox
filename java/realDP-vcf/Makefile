SHELL=/bin/bash
HOME=$(shell echo $$HOME)
PWD=$(shell pwd)
BAMS=../../bam.examples
PICARD_JARS=$(PICARD)
SAM_JAR=$(shell ls $(PICARD_JARS)/*sam-*.jar)
PICARD_JAR=$(shell ls $(PICARD_JARS)/*picard-*.jar)
CP=-classpath "$(SAM_JAR):$(PICARD_JAR):."
TOOL_NAME=vcfAddCoverage

all: output.vcf

output.vcf: $(TOOL_NAME).class
	cat merged.vcf | java -Xmx1g $(CP) $(<:.class=)

$(TOOL_NAME).class: $(TOOL_NAME).java
	javac $(CP) $<
