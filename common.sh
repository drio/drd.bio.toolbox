#!/bin/bash
#
set -e
#set -x

main_dir="/data/next_gen_1/drio_scratch/macaque"
bin_dir="/data/next_gen_1/drio_scratch/macaque/bin"
dnaa="$bin_dir/dnaa"
dwgsim="$dnaa/dwgsim/dwgsim"
novo_bin_dir="$HOME/projects/novoalign/bin/macosx"
novo_bf_2_novo="`dirname ${BASH_SOURCE[0]}`/third-party/bfast2novo.pl"
true_snps="true.snps.txt"
novo_reads="novo.read1.fastq novo.read2.fastq"
e_rate="0.05"
bwa_read1="bwa.read1.fastq"
bwa_read2="bwa.read2.fastq"
bf_reads="bf.reads"

novo_index_ss="$bin_dir/novoalign/bin/linux/novocraft/novoindex"
novo_index_cs="$bin_dir/novoalign/bin/linux//novoalignCS/novoindex"
novo_ss_index="novo.index.ss"
novo_cs_index="novo.index.cs"

bwa="$bin_dir/bwa/bwa"
