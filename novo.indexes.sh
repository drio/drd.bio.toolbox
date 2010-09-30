#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

ff=$1

cat <<-EOF
#
# Usage:
# <script> <fast_file_ref_genome>
#
# NOTE: Novoalign will spawn as many threads as core(s) available + 1 !!
#
$novo_index_ss -k 15 -s 2 $novo_ss_index.`basename $1` $ff
$novo_index_cs -k 15 -s 2 -c $novo_cs_index.`basename $1` $ff
EOF
