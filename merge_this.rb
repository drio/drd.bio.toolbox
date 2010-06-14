#!/usr/bin/env ruby19
#
# Merges bams
#
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

puts DATA.read.gsub(/INPUTS/, n_i).gsub(/XXOUTXX/, out_fn)

__END__
bsub -o out -e err \
-q high \
-R 'rusage[mem=8000]span[hosts=1]' \
-J merge_this.XXOUTXX \
/stornext/snfs1/next-gen/software/jdk1.6.0_01/bin/java -jar -Xmx8g \
/stornext/snfs1/next-gen/software/picard-tools/current/MergeSamFiles.jar \
INPUTS \
TMP_DIR=/space1/tmp \
OUTPUT=XXOUTXX \
VALIDATION_STRINGENCY=SILENT 
#USE_THREADING=true \
#VALIDATION_STRINGENCY=STRICT
