SHELL=/bin/bash
HOME=$(shell echo $$HOME)
PWD=$(shell pwd)
BAMS=../../bam.examples
PICARD_JARS=$(PICARD)
SAM_JAR=$(shell ls $(PICARD_JARS)/*sam-*.jar)
PICARD_JAR=$(shell ls $(PICARD_JARS)/*picard-*.jar)
CP=-classpath "$(SAM_JAR):$(PICARD_JAR):$(PWD)"
FLAGS=-Xlint
TOOL_NAME=vcfAddCoverage

all: test

test: test/merged.vcf RCoverage.class RealRDP.class
	cd test; java -Xmx1g $(CP) -ea RealRDP merged.vcf

#output.vcf: $(TOOL_NAME).class
#	cat data/merged.vcf | java -Xmx1g $(CP) $(<:.class=)

%.class: %.java
	javac $(FLAGS) $(CP) $<
