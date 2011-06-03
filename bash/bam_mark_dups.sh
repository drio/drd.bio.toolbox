#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

usage()
{
  [ ".$1" != "." ] && echo "ERROR: $1"
  echo "Usage:"
  echo "$0 <output file> <input sam|bam>"
  exit 0
}

run()
{
  j=${root}.$2
  cmd="$1"
  echo $1
}

of=$1
if=$2

[ ".$if" == "." ]  && usage "I need a sam file"
[ ".$of" == "." ] && usage "I need an output file"
[ ! -d "$picard_jars_dir" ]  && usage "Couldn't find picard dir: $picard_jars_dir"

#bi="$root.bam"
#bs="$root.sort.bam"
#bd="$root.sort.dups.bam"

#if [ `echo ${if#*.} | grep "sam$"` ]
#then
#  echo "samtools view -hbS $ref $if > $bi"
#else
#  bi=$if
#fi

cat <<EOF
$java -jar -Xmx26g $picard_jars_dir/MarkDuplicates.jar \
TMP_DIR=$tmp \
METRICS_FILE='metric_file.picard' \
INPUT=$if \
OUTPUT=$of \
MAX_RECORDS_IN_RAM=2500000 \
VALIDATION_STRINGENCY=$PICARD_VALIDATION
EOF
