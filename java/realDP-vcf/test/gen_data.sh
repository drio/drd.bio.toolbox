#!/bin/bash
#
set -e

s_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
b_dir="/stornext/snfs6/rogers/drio_scratch/rhesus/crv.v2/wes"
#ids="34601 32510 40459"
#ids="34601 32510 40459 34608 35117 35058 32738 35048 31506 32509 23138 34611 36070 34604 34597 34612"
ids="34601 32510 40459 34608 35117 35058 32738 35048 31506 32509 23138 34611 36070 34604 34597 34612 34607 34594 27347 34051"
n_lines=100000
REF="/stornext/snfs6/rogers/drio_scratch/genomes/rhemac2.indian_macaque_no_phix.fixed.fa"

cd $s_dir
echo ">> Generating test datasets."
for i in $ids
do
  bam=${i}.bam
  rm -f $bam $i.vcf* $bam.bai
  samtools view -h "$b_dir/${i}.bam" | head -${n_lines} | samtools view -Sb - > $bam
  samtools index $bam

  o=$i.vcf.gz
  samtools mpileup -R -E -uf $REF $bam | bcftools view -vcg - | bgzip -c > $o
  tabix -p vcf $o
  vcfs="$vcfs $o"
done
vcf-merge $vcfs > merged.vcf

