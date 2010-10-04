#!/bin/bash
#

cat<<EOF
BWA
  $ solid2fastq.pl 0044_20100421_2_SP_ANG_AUT_N_200716041_1_1sA_01003310991_3_ bwa
  $ gzip -d bwa.read1.fastq.gz &
  $ gzip -d bwa.read2.fastq.gz &
  $ wait
  $ mkdir -p splits.read1
  $ mkdir -p splits.read2
  $ split -a 1 -d -l 40000000 ../bwa.read1.fastq ./splits.read1/read1.fastq. &
  $ split -a 1 -d -l 40000000 ../bwa.read2.fastq ./splits.read2/read2.fastq. &
  $ wait
  $ bwa aln -c -t8 ref fq > i.sai
  $ bwa samse ref i.sai fq > i.sam
  $ samtools view -bS i.sam > i.bam
  $ samtools flagstat i.bam  > i.stats.txt

NOVOALIGN
  $ novo_ss_dir/novoindex novo_index_ss ff > ./novo.index.ss.txt
  $ novo_cs_dir/novoindex -c novo_index_cs ff > ./novo.index.cs.txt
  $ novo_ss_dir/novoalign   -c threads -d novo_index_ss -f novo_reads -o SAM > novo.sam 2> ./novo.align.ss.txt
  $ novo_cs_dir/novoalignCS -c threads -d novo_index_cs -f novo_reads -o SAM > novo.sam 2> ./novo.align.cs.txt

BFAST
  $ bfast_bin  match -T"/tmp/" -n threads -A1 -l -f "./bfast.indexes/ff" -r bf_reads > bf.matches.bmf 2>/dev/null
  $ bfast_bin  localalign -n threads -A1  -f prefix -m bf.matches.bmf > bf.aligned.baf 2>/dev/null
  $ # postprocess
  $ # -O1: output SAM; 
  $ # -a3: Choose uniquely the alignment with the best score
  $ # -U : Don't perform pairing
  $ bfast_bin postprocess -f "./bfast.indexes/ff" -U -i bf.aligned.baf -n threads -a3 -O1 > bf.sam 2>/dev/null

-------------------------------------
DIST dir: `dirname ${BASH_SOURCE[0]}`
-------------------------------------
EOF

ls -l `dirname ${BASH_SOURCE[0]}`
