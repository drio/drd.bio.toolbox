#!/bin/bash
#
set -e 

usage()
{
  cat<<-EOF
Usage:
  script <bwa_dwgsim_eval> <bf_dwgsim_eval> <novo_dwgsim_eval> <plot title> <plot_output_name>

Examples:
$ run_R_roc_accuracy.sh all.accuracy.txt roc.accuracy.pdf
$ for i in snps inss dels; do w="\$i" && run_R_roc_snps.sh bwa*.\$w bf*.\$w novo*.\$w  SNPs roc.\$w.pdf; done  
EOF
}

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

d_path=`dirname ${BASH_SOURCE[0]}`
r_script="${d_path}/../R/roc.snps.R"

[ ! -f $1 ] && usage && error "bwa dwgsim eval file not found"
[ ! -f $2 ] && usage && error "bf dwgsim eval file not found"
[ ! -f $3 ] && usage && error "novo dwgsim eval file not found"
[ ".$4" == "." ] && usage && error "Title not provided"
[ ".$5" == "." ] && usage && error "Output name not provided"

R CMD BATCH "--vanilla --args in.bwa='$1' in.bf='$2' in.novo='$3' in.o_file='$5' in.title='$4'" $r_script
echo -e "bwa:\t$1"
echo -e "bf:\t$2"
echo -e "novo:\t$3"
ls -lach $5
