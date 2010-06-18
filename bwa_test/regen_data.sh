#!/bin/bash

rm -rf raw_v4 raw

mkdir -p raw
cd raw
for i in /stornext/snfs1/next-gen/drio-scratch/bfast_related/bwaaln.testing/ANG_RC_LCA_K86_3_1sA_01003280864_2/raw/*
do 
  head -200001 $i > `basename $i` 
done
cd ..

mkdir -p raw_v4
cd raw_v4
for i in /stornext/snfs1/next-gen/drio-scratch/bfast_related/bwaaln.testing/ANG_AUT_N_200716041_1_1sA_01003310991_3/raw/*
do 
  head -200001 $i > `basename $i` 
done
