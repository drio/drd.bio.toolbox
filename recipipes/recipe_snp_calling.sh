#!/bin/bash
set -e 
#set -x 

source "`dirname ${BASH_SOURCE[0]}`/../bash/common.sh"

log()
{
  echo ""
  echo "### `date` >> $1"
}

usage()
{
  [ ".$1" != "." ] && echo "ERROR: $1"
  echo "usage: `basename $0` -f fasta_ref_file <bam1 .. bamN>"
  exit 0
}

while getopts "f:h" OPTION
do
  case $OPTION in
    f)
      fasta_file=$OPTARG
      ;;
    ?)
      usage
      exit 2
      ;;
  esac
done
shift $(($OPTIND - 1))

[ ".$fasta_file" = "." ] && usage "Need a fasta file."
[ ! -f "$fasta_file" ] && usage "Fasta file: <$fasta_file> not found"
input_bams=""
for b in $@; do
  [ ! -f $b ] && usage "bam: <$b> not found."
  input_bams="$input_bams $b" 
done


log "Merge bams ($input_bams)"
merged_bam="merged.bam"
merge_this $input_bams $merged_bam | bash

log "mark dups"
merged_dups_bam="merged.dups.bam"
bam_mark_dups.sh $merged_dups_bam $merged_bam | bash

base_cov="base_coverage.txt"
log "calculating base coverage"
std_pileup $merged_dups_bam | awk '{if($3>1) print;}' > $base_cov

log "pileup"
pileup.this.sh $merged_dups_bam $fasta_file | bash

log "dumping stats"
output="snp_stats.json"
snp_stats.rb $base_cov *.var_filter.qual_threshold > $output

log "Done."
