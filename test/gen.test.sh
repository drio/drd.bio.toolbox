#!/bin/bash
#
# Careful with the testing
# Bwa asigns a 0 to reads that are ambiguous
# Bfast tosses those ones (new realease has a behaviour that mimics bwa)

#b="/stornext/snfs1/next-gen/drio-scratch/bfast_related/bwaaln.testing/ANG_AUT_N_200716041_1_1sA_01003310991_3/bf2/50_25_n3/merged/01003310991_3.bf2.n3.merged.bam"
b="/stornext/snfs1/next-gen/drio-scratch/bfast_related/bwaaln.testing/ANG_AUT_N_200716041_1_1sA_01003310991_3/bwa/n3/sampe/01003310991_3.bwa.n3.sampe.bam"

samtools view -h $b | head -1000000 | samtools view -h -S - > i.sam
