#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

usage()
{
  [ ".$1" != "." ] && echo "ERROR: $1"
  echo "Usage:"
  echo "$0 <input_bam> <output_seed> <fasta_ref>"
  exit 0
}

bf=$1
seed=$2
ref=$3

[ ".$bf" == "." ]  && usage "I need a bam"
[ ".$seed" == "." ] && usage "I need the output seed"
[ ".$ref" == "." ] && usage "I need a ref genome (fasta)"
[ ! -f "$bf" ]  && usage "Can't find bam."
[ ! -f "$ref" ]  && usage "Can't find ref."

cat <<-EOF
# -m INT      minimum gapped reads for indel candidates [1] 
# -F FLOAT    minimum fraction of gapped reads for candidates [0.002]
# -D          output per-sample DP
# -S          output per-sample SP (strand bias P-value, slow)
# -g          generate BCF output 
# -f FILE     reference sequence file
# -u          uncompress output
#
# -b        output BCF instead of VCF 
# -v        output potential variant sites only (force -c)
# -c        SNP calling
# -g        call genotypes at variant sites (force -c)
#
samtools mpileup -uf $ref $bf | bcftools view -bvcg - > var.raw.${seed}.bcf
# -d INT    minimum read depth
# -D INT    maximum read depth
bcftools view var.raw.${seed}.bcf | vcfutils.pl varFilter -D100 -d3 > var.filtered.${seed}.vcf
EOF
