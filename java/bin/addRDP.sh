#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

[ -t 0 ] && echo "Need data in stdin." && exit 1

echo $DIR
cat - | \
java -Xmx4g \
-classpath "$PICARD/sam-1.79.jar:$PICARD/picard-1.79.jar:$DIR/../realDP-vcf" vcfAddCoverage
