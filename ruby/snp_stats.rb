#!/bin/env ruby
#
require File.expand_path(File.dirname(__FILE__)) + "/common"

include Help
usage =<<EOF
Given the base coverage and varfiltered samtools output
generates basic stats.

Usage: 
  $ #{File.basename(__FILE__)}" <base_coverage.txt> <var_filtered_output>
EOF
Help.set_usage_text usage

Help.error "Invalid number of arguments." if ARGV.size != 2
cov_fn, pu_fn = ARGV
Help.error "coverage file not found. " unless File.exists?(cov_fn)
Help.error "var_filtered pu file not found. " unless File.exists?(pu_fn)

# callable positions (we have at least two reads covering the locus)
callable = 0
Common::open_file(cov_fn).each_line do |l|
  callable += 1 if l.split[2].to_i > 1
end

# Find number of snps and indels (from var_filtered pileup)
subs = indels = 0
File.open(pu_fn).each_line do |l|
  s = l.split
  if s[2] == '*'
    indels += 1
  else
    subs += 1
  end
end

# Dump results
total_events = subs + indels
puts <<EOF
{
  callable: #{callable},
  substitutions: #{subs},
  indels: #{indels},
  total_events: #{total_events},
  percentage_difference: #{(total_events.to_f * 100) / callable.to_f},
}
EOF
