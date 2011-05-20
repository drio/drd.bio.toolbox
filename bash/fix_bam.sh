#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

usage()
{
  [ ".$1" != "." ] && echo "ERROR: $1"
  echo "Usage:"
  echo "$0 <bam>"
  exit 0
}

if=$1

[ ".$if" == "." ] && usage "I need an input_bam"
[ ! -d "$ipipe_java" ]  && usage "Couldn't find illumina pipe dir: $ipipe_java"

cat <<EOF
$java -jar -Xmx4g $ipipe_java/FixCIGAR.jar \
I=$if \
TMP_DIR=$tmp \
VERBOSITY=ERROR \
VALIDATION_STRINGENCY=STRICT
EOF
