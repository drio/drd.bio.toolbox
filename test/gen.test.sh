#!/bin/bash
#

samtools view -h /stornext/snfs4/next-gen/solid/analysis/solid0513/2010/04/0513_20100421_1_SP_ANG_AUT_N_16685085_1_1sA_01003310985_4/output/0513_20100421_1_SP_ANG_AUT_N_16685085_1_1sA_01003310985_4.sorted.dups.bam | head -100000 | samtools view  -h -S - > i.sam
