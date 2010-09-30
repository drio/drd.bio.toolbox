#!/bin/bash
#
set -e
#set -x

root=$1
if=$2
java="/stornext/snfs1/next-gen/software/jdk1.6.0_01/bin/java"
picard="/stornext/snfs1/next-gen/software/picard-tools/current"
tmp="/space1/tmp"

usage()
{
  [ ".$1" != "." ] && echo "ERROR: $1"
  echo "Usage:"
  echo "$0 <root seed> <sam|bam file>"
  exit 0
}

run()
{
  j=${root}.$2
  cmd="$1"
  #echo "bsub -q high -R 'rusage[mem=8000]' -J $j \"$cmd\""
  echo $1
}

[ ".$if"  == "." ]  && usage "I need a sam file"
[ ".$root" == "." ] && usage "I need a root seed"
[ ! -d "$picard" ]  && usage "Couldn't find picard dir: $picard"

bi="$root.bam"
bs="$root.sort.bam"
bd="$root.sort.dups.bam"

cat <<EOF
#!/bin/bash
#
#BSUB -J $root.bam.sort.dup 
#BSUB -o lsf.out
#BSUB -e lsf.err
#BSUB -q normal
#BSUB -R 'rusage[mem=8000]'
#BSUB -u "deiros@bcm.edu"
#
set -e
EOF

if [ `echo ${if#*.} | grep "sam$"` ]
then
  echo "samtools view -hbS $ref $if > $bi"
else
  bi=$if
fi

# 1. Sort by Coordinate
cat <<EOF
$java -jar -Xmx2g $picard/SortSam.jar \
TMP_DIR=$tmp \
INPUT=$bi \
OUTPUT=$bs \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=STRICT
#VALIDATION_STRINGENCY=SILENT
EOF

# 2. Mark dups
cat <<EOF
$java -jar -Xmx6g $picard/MarkDuplicates.jar \
TMP_DIR=$tmp \
INPUT=$bs \
OUTPUT=$bd \
METRICS_FILE='/tmp/metric_file.picard' \
VERBOSITY=ERROR \
VALIDATION_STRINGENCY=STRICT
#VALIDATION_STRINGENCY=SILENT
EOF
