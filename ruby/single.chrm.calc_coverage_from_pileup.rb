#!/bin/env ruby
#
# Finds coverage using samtools pileup input.
#
require File.expand_path(File.dirname(__FILE__)) + "/common"
#require "profile"

include Help
include Common
include Statistics
usage=<<-EOF
Usage: #{File.basename(__FILE__)} <genome_or_scaffold_or_contig size>

Example: 

  $ zcat pu.chrm10.gz | awk '{print $1" "$2" "$8}' | uniq | \
  calc_coverage_stats 160000000
  # size of Gibbon chrm10 in assembly (Not actual value)

Input:
chrm position read_coverage

Note   : The tools expects the bam to be sorted by coordinate.
WARNING: Also, we use uniq to skip join indels!
EOF
Help.set_usage_text usage

# Same as process_pileup, but does not stores anything
# in the array, it dumps directly.
module PileUpProcessing
  def self.process_pu_and_dump(gsize)
    gp     = 0   # global position
    offset = 0   # pileup uses positions relative to the scaffolds.
    $stdin.each_line do |l_stdin|
      s, p, c = l_stdin.split # scaffold name, position, coverage
      p = p.to_i
      while gp < p # We don't have coverage for that region
        printf "%s %s\n", gp, 0; gp += 1
      end
      printf "%s %s\n", gp, c # We have coverage in pileup, use it
      gp += 1
    end
    # If there are more locations without coverage, dump them
    while gp < gsize 
      printf "%s %s\n", gp, 0; gp += 1
    end
  end
end

# Main
# 
if __FILE__ == $0
  Help.error "Invalid number of arguments." if ARGV.size != 1
  Help.error "No data in stdin." unless Common.data_in_stdin?
  g_size = ARGV[0].to_i
  PileUpProcessing::process_pu_and_dump g_size
end
