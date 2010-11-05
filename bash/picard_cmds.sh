#!/bin/bash
#

source "`dirname ${BASH_SOURCE[0]}`/common.sh"
j="java -Xmx4g -jar"
tmp="/tmp"

cat <<-EOF
#
# Usage:
# <script> <bam_file> <reference_file>
#
EOF

echo "set -e"
r="QualityScoreDistribution"
echo $j $PICARD/$r.jar INPUT=$1 OUTPUT=$r.txt CHART_OUTPUT=$r.pdf TMP_DIR=$tmp 

r="CollectAlignmentSummaryMetrics"
echo $j $PICARD/$r.jar INPUT=$1 OUTPUT=$r.txt TMP_DIR=$tmp REFERENCE_SEQUENCE=$2

r="CollectGcBiasMetrics"
echo $j $PICARD/$r.jar INPUT=$1 OUTPUT=$r.txt CHART_OUTPUT=$r.pdf TMP_DIR=$tmp REFERENCE_SEQUENCE=$2

r="CollectInsertSizeMetrics"
echo $j $PICARD/$r.jar INPUT=$1 OUTPUT=$r.txt  TMP_DIR=$tmp HISTOGRAM_FILE=$r.pdf

r="MeanQualityByCycle"
echo $j $PICARD/$r.jar INPUT=$1 OUTPUT=$r.txt CHART_OUTPUT=$r.pdf TMP_DIR=$tmp

r="QualityScoreDistribution"
echo $j $PICARD/$r.jar INPUT=$1 OUTPUT=$r.txt CHART_OUTPUT=$r.pdf TMP_DIR=$tmp
