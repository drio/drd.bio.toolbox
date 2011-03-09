#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

root=$1
if=$2

usage()
{
  [ ".$1" != "." ] && echo "ERROR: $1"
  echo "Usage:"
  echo "$0 <root seed output file> <sam|bam file>"
  exit 0
}

run()
{
  j=${root}.$2
  cmd="$1"
  echo $1
}

[ ".$if"  == "." ]  && usage "I need a sam file"
[ ".$root" == "." ] && usage "I need a root seed"
[ ! -d "$picard_jars_dir" ]  && usage "Couldn't find picard dir: $picard_jars_dir"

bi="$root.bam"
bs="$root.sort.bam"
bd="$root.sort.dups.bam"

if [ `echo ${if#*.} | grep "sam$"` ]
then
  echo "samtools view -hbS $ref $if > $bi"
else
  bi=$if
fi

cat <<EOF
$java -jar -Xmx6g $picard_jars_dir/MarkDuplicates.jar \
TMP_DIR=$tmp \
INPUT=$if \
OUTPUT=$bd \
METRICS_FILE='/tmp/metric_file.picard' \
VERBOSITY=ERROR \
VALIDATION_STRINGENCY=STRICT
#VALIDATION_STRINGENCY=SILENT
EOF
