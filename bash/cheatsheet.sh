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

Pipeline 1: Generation of ROC curves to test performance of aligners with simulated data (CS)

  0.  bash/*.indexes.sh              : Generate indexes for your aligner of choice.
  1.  bash/calculate_coverage.sh     : Find out how many reads you need to simulated your coverage.
  2.  bash/simulate_reads.sh         : Simulated reads for different aligners
  3.  bash/split_reads.sh            : Split reads 
  4.  bash/align_with_*.sh           : Align with your aligner of choice.
  5.  ruby/merge_this.rb             : Merge the generated SAMs (generate BAM)
  6.  bash/bam_sort_dups.sh          : Sort and mark dups. 
  7.  dwgsim_eval                    : Evaluate the quality of the alignments.
  8.  ruby/gen_accuracy_plot_data.rb : Prepare dwgsim_eval data for plotting.
  9.  bash/pileup.this.sh            : Generate pileup from BAM.
  10. dwgsim_pileup_eval.pl          : Evaluate the snp calls (pileup). 
  11. bash/gen_eval_data_for_plot.sh : Prepare dwgsim_eval data output for plotting.
  12. bash/run_R_roc_accuracy.sh     : Generate ROC curve for accuracy.
  13. bash/run_R_roc_snps.sh         : Generate ROC curves for SNP/indel detection.

EOF

echo "---------"
echo "You may want to add this to your .bashrc (if you use bash): "
echo "source `dirname ${BASH_SOURCE[0]}`/common.sh"
echo 'export PATH=$PATH:$dnaa/dwgsim:$novo_cs:$novo_ss:$bin_dir/bwa:$bin_dir/bfast/bfast'

