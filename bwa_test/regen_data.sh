#!/bin/bash

rm -rf raw*
rm -f *.sai *.txt *.gz *.sam *.bam *.fastq *.png

mkdir -p raw
cd raw
for i in /stornext/snfs1/next-gen/drio-scratch/bfast_related/bwaaln.testing/ANG_RC_LCA_K86_3_1sA_01003280864_2/raw/*
do 
  head -100001 $i > `basename $i` 
done
cd ..

mkdir -p raw_v4
cd raw_v4
for i in /stornext/snfs1/next-gen/drio-scratch/bfast_related/bwaaln.testing/ANG_AUT_N_200716041_1_1sA_01003310991_3/raw/*
do 
  head -100001 $i > `basename $i` 
done
cd ..

cp -r raw_v4 raw_v4_fixed
cd raw_v4_fixed
mv 0044_20100421_2_SP_ANG_AUT_N_200716041_1_1sA_01003310991_3_F5-P2.csfasta 0044_20100421_2_SP_ANG_AUT_N_200716041_1_1sA_01003310991_3_R3.csfasta
mv 0044_20100421_2_SP_ANG_AUT_N_200716041_1_1sA_01003310991_3_F5-P2_QV.qual 0044_20100421_2_SP_ANG_AUT_N_200716041_1_1sA_01003310991_3_R3_QV.qual
cd ..
