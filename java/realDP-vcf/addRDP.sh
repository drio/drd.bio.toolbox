#!/bin/bash
#
set -e
#set -x

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export CLASSPATH="$DIR/bin:$DIR/picard/picard-1.56.jar:$DIR/picard/sam-1.56.jar:."
java AddRDP $@
