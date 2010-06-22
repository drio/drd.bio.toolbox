#!/bin/bash

r="0044_20100421_2_SP_ANG_AUT_N_200716041_1_1sA_01003310991_3_"
p="./raw_v4_fixed"
ref="/stornext/snfs4/next-gen/solid/bwa.references/h/hsap.36.1.hg18/hsap_36.1_hg18.cs"

run_bwa.rb -r $p -u $r -f $ref
