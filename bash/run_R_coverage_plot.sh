#!/bin/bash
#
#set -e 

#script <output.png> <window_size> <title> <pileup_file>
usage()
{
  cat<<-EOF
Usage:
  script <output.png> <title>

Example:
$ cat i.txt | awk '{print $2" "$8 }' | uniq |awk '{print $2}'| calc_coverage_stats 50 | run_R_coverage_plot.sh cov.png title_here
$ samtools pileup my.bam | awk '{print $2" "$8 }' | uniq |awk '{print $2}'| calc_coverage_stats 50 | run_R_coverage_plot.sh cov.png title_here
EOF
}

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

d_path=`dirname ${BASH_SOURCE[0]}`
#r_script="${d_path}/../R/coverage.plot.R"
r_script="${d_path}/../R/coverage.plot.mean.version.R"

#if [ ".$1" == "." ] || [ ".$2" == "." ] || [ ".$3" == "." ] || [ ! -f $4 ]
if [ ".$1" == "." ] || [ ".$2" == "." ]
then
  usage
  exit 1
fi

#awk '{print $2" "$8}' $4 | R CMD BATCH "--vanilla --slave --args in.o_file='$1' in.window='$2' in.title='$3'" $r_script
R CMD BATCH "--vanilla --slave --args in.o_file='$1' in.title='$2'" $r_script
