#!/bin/bash
#
set -e

#bw_path="/tmp/drd/breakway.0.7"
#bw_path="/data/rogers/drio_scratch/bb/local/breakway.0.7"
#bam="/data/rogers2/drio_scratch/gibbon/AWG_NLEU.00_000pA/AWG_NLEU.00_000pA.merged.nodups.q_gt_30.bam"
#ff="/data/rogers2/drio_scratch/genomes/Nleu1.0.fixed.fa"
#fai="/data/rogers2/drio_scratch/genomes/Nleu1.0.fixed.fa.fai"
#rl=50 # read length
dt_output="./dtrans.txt"
expect_window=4000
dt_inc=100
dt_bin=500
few_reads="100000" # when testing, use this amount of reads
#
# minimum number of reads in a bin in the Dtranslocations file for it to be included in analysis
cut_off=4
# is the minimum score threshold that events must exceed in order to be included in the output.
score=0
# --strand: can have one of two values: "SAME" or "OPPOSITE". These refer to the strand
# orientation of reads in your data set. If a normal paired read from your data
# set has both ends of the read on the same strand (+/+ or -/-), use "SAME", and
# if they're on opposite strands (+/- or -/+), use "OPPOSITE".
# This is determined by your sequencing technology and aligner. For reference,
# Illumina GAII data is typically "OPPOSITE", while ABI SOLiD data is typically
# "SAME". If you are not sure which your data is, you can extract a subset of the
# aligned reads and determine which strand orientation dominates your data set to
# find out.
strand=""
# --ordering: can have one of two values: "12" or "21". The digits refer to the
# read number, as paired reads are always paired up as "read 1" and "read 2".
# If read 2 is downstream of read 1 (5'-read1-read2-3'), --ordering is 12. If
# read 1 is downstream of read 2 (5'-read2-read1-3'), --ordering is 21. For
# reference, Illumina GAII data is typically "12", while ABI SOLiD data is
# typically "21". Again, you can determine from your own dataset what the
# ordering is by extracting a subset of aligned reads and seeing which ordering
# dominates.
ordering=""
testing=0 # By default, no testing mode.

usage()
{
  [ ".$1" != "." ] && echo "ERROR: $1"
  cat <<-EOF
usage: `basename $0` options

  -p <path>  : path to breakway
  -b <bam>   : full path do bam file.
  -f <fasta> : full path do ref file.
  -r <rl>    : read length.
  -s <strand>: SAME or OPPOSITE (see: breakway compendium).
               SAME (SOLiD), OPPOSITE: Illumina
  -o <order> : 21 12 (see: breakway compendium)
               12: Illumina
               21: SOLiD
  [-t]       : Enable testing mode
  [-g]       : Genome size

Example(s):
  $ `basename $0` -p ./breakway.0.7 -b i.bam -f fasta.fa -r 50 -s SAME -o 21 -t
EOF
  exit 0
}

log()
{
  echo "`date` >> $1" >&2
}

testing_reads()
{
  log "WARNING (test mode enabled); extracting few reads ($few_reads)"
  samtools view -h $bam | head -${few_reads} | samtools view -bS - > i.bam 2>/dev/null
  bam="i.bam"
}

find_ped_stats()
{
  log "finding PED stats [./pen.txt]"
  perl $bw_path/scripts/breakway.parameters.pl 0 $expect_window $bam > pen.txt 2> /dev/null
  ped_mean=`cat pen.txt  | grep Mean     | awk  -F= '{print $2}' | sed 's/ //g'`
  ped_std=`cat pen.txt   | grep Standard | awk  -F= '{print $2}' | sed 's/ //g'`
  ped_lower=`cat pen.txt | grep lower    | awk  -F= '{print $2}' | sed 's/ //g'`
  ped_upper=`cat pen.txt | grep upper    | awk  -F= '{print $2}' | sed 's/ //g'`
  ped_lower=$[$ped_lower+100]
  ped_upper=$[$ped_upper+100]
  log "PED mean     : ${ped_mean}"
  log "PED std      : ${ped_std}"
  log "PED 95% upper: ${ped_upper}"
  log "PED 95% lower: ${ped_lower}"
}

get_bam_stats()
{
  log "bam stats [./stats.txt]"
  dbamstats $bam > stats.txt 2>/dev/null
  nm_pe_reads=`cat stats.txt | grep "mapped paired reads:" | awk -F: '{print $2}'`
  min_cs=`echo "$ped_lower/$dt_inc" | bc`
  max_cs=`echo "$ped_upper/$dt_inc" | bc`
  hap_ccov=`echo "scale=2;(($nm_pe_reads*(2*$rl+$ped_mean)) / $g_size) / 2" | bc`
  std_ccov=`echo "scale=2;$nm_pe_reads*$ped_std/$g_size" | bc`
  log "# mapped pe reads            : ${nm_pe_reads}"
  log "haploid clone coverage       : ${hap_ccov}x"
  log "st dev haploid clone coverage: ${std_ccov}"
  log "min cluster size             : ${min_cs}"
  log "max cluster size             : ${max_cs}"
}

run_dtranslocations()
{
  log "running dtranslocations [${dt_output}]"
  dtranslocations -i $dt_inc -b $dt_bin -l $ped_lower -L $ped_upper $bam > $dt_output 2>/dev/null
}

find_genome_size()
{
  log "finding genome size"
  if [ ".$g_size" == "." ]; then
    log "Calculating genome size"
    g_size=`cat $ff | ruby -ne 'BEGIN{@t=0}; @t = @t + $_.chomp.size unless $_ =~ /^>/; END{puts @t.to_i}'`
  else 
    log "User provided genome size."
  fi
  log "genome size: $g_size bp"
}

find_svs()
{
  log "finding sv [./breakway.out.bw]"
  #cat <<-EOF
  perl $bw_path/breakway.run.pl \
  --dtransfile $dt_output \
  --mincs $min_cs --maxcs $max_cs \
  --interval $dt_inc \
  --cutoff $cut_off \
  --bamfile $bam \
  --fai $fai \
  --filehandle breakway \
  --bwfolder $bw_path \
  --mindist $ped_lower --maxdist $ped_upper --mean $hap_ccov --stdev $std_ccov --score $score \
  --strand $strand --ordering $ordering \
  > breakway.out.bw 2> breakway.stderr
  #EOF
  #--segdupsfile [genomicSuperDups.txt] --selfchainfile [chainSelf.txt] --repmaskerfile [rmsk.txt]
}

#
# Main
while getopts "p:b:f:r:s:o:g:t" OPTION
do
  case $OPTION in
    p)
      bw_path=$OPTARG
      ;;
    b)
      bam=$OPTARG
      ;;
    f)
      ff=$OPTARG
      ;;
    r)
      rl=$OPTARG
      ;;
    s)
      strand=$OPTARG
      ;;
    o)
      ordering=$OPTARG
      ;;
    g)
      g_size=$OPTARG
      ;;
    t)
      testing=1
      ;;
    h)
      usage
      exit 0
      ;;
    ?)
      usage
      exit 2
      ;;
  esac
done

fai="${ff}.fai"
bai="${bam}.bai"

# Double check input
#
[ ! -d $bw_path   ]   && usage "Cannot find breakway dir."
[ ! -f $bam       ]   && usage "Cannot find bam file."
[ ! -f $bai       ]   && usage "Bam index not found: $bai"
[ ! -f $ff        ]   && usage "ff not found."
[ ! -f $fai       ]   && usage "Index of reference not found: $fai"
[ "$rl" == ".$rl" ]   && usage "Need to know the read length."
[[ $strand != "SAME" && $strand != "OPPOSITE" ]] && usage "Incorrect strand." 
[[ $ordering != "12" && $ordering != "21" ]] && usage "Incorrect ordering." 

[ $testing == "1" ] && testing_reads

find_genome_size
find_ped_stats
get_bam_stats
run_dtranslocations
find_svs
log "done"
