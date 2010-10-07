#!/bin/env ruby
#
# Merges bams
#
require File.expand_path(File.dirname(__FILE__)) + "/common"

def usage
  puts <<END
Usage: #{__FILE__} <list of files to merge> <out_bam>
END
  exit 1
end


ARGV.size < 2 && usage
files_to_merge = ARGV[0..-2]
out_fn = ARGV[-1]

n_i = files_to_merge.inject("") do |f, i| 
  f << "INPUT=" + i + " "
end

cmd = DATA.read.gsub(/_INPUTS_/, n_i).gsub(/_OUT_/, out_fn)
ch = Common::load_common
cmd.gsub!(/_JAVA_/, ch["java"])
cmd.gsub!(/_JAR_/, ch["picard_jars_dir"] + "/MergeSamFiles.jar")
cmd.gsub!(/_TMP_/, ch["tmp"])
puts cmd
#`#{cmd}`

__END__
#bsub -o out -e err \
#-q high \
#-R 'rusage[mem=8000]span[hosts=1]' \
#-J merge_this.XXOUTXX \
#/stornext/snfs1/next-gen/software/picard-tools/current/MergeSamFiles.jar \
_JAVA_ -jar -Xmx8g _JAR_ \
_INPUTS_ \
TMP_DIR=_TMP_ \
OUTPUT=_OUT_ \
VALIDATION_STRINGENCY=SILENT 
#USE_THREADING=true \
#VALIDATION_STRINGENCY=STRICT
