#!/bin/bash
#
#set -e
#set -x

error()
{
  echo "ERROR: $1"
  exit 1
}

script_this()
{
(
  cat<<EOF
#!/bin/bash
cd `pwd`
$1
EOF
) > $2
}

main_dir="/data/next_gen_1/drio_scratch/macaque"
bin_dir="/data/next_gen_1/drio_scratch/macaque/bin"
dnaa="$bin_dir/dnaa"
dwgsim="$dnaa/dwgsim/dwgsim"
novo_bin_dir="$HOME/projects/novoalign/bin/linux"
novo_bf_2_novo="`dirname ${BASH_SOURCE[0]}`/../third-party/bfast2novo.drd.pl"
true_snps="true.snps.txt"
novo_reads="novo.read1.fastq novo.read2.fastq"
e_rate="0.05"
bwa_read1="bwa.read1.fastq"
bwa_read2="bwa.read2.fastq"
bf_reads="bf.reads"
threads="4"

novo_index_ss="$bin_dir/novoalign/bin/linux/novocraft/novoindex"
novo_index_cs="$bin_dir/novoalign/bin/linux//novoalignCS/novoindex"
novo_ss_index="novo.index.ss"
novo_cs_index="novo.index.cs"
novo_cs="$bin_dir/novoalign/bin/linux/novoalignCS"
novo_ss="$bin_dir/novoalign/bin/linux/novocraft"


bwa="$bin_dir/bwa/bwa"
tmp="/space1/tmp/"

java="java"
R="R"
picard_jars_dir="/data/next_gen_1/drio_scratch/software/picard-tools"
ruby="/data/next_gen_1/drio_scratch/software/bin/ruby"

samtools="samtools"
samtools_perl="samtools.pl"

submit_bin="msub"  # MOAB
#submit_bin="qsub" # PBS

# Bfast
#
bfast_bin="$bin_dir/bfast/bfast/bfast"
h_masks=(
1111111111111111111111
111110100111110011111111111
10111111011001100011111000111111
11111111100101111000001100011111011
1111111110001111110011111111
111111011010011000011000110011111111
11111111111110011101111111
1111011000011111111001111011111
11110110001011010011100101111101111
1111111001000110001011100110001100011111
)

