#!/bin/env ruby
#
# Converts flagstats output to csv
# 
#
require File.expand_path(File.dirname(__FILE__)) + "/common"
include Common
include Help
usage =<<EOF
Usage: 
  $ dbamstats my.bam | #{File.basename(__FILE__)}"
EOF
Help.set_usage_text usage
Help.error "No data in stdin." unless data_in_stdin?

#Number of raw paired reads:88004	
#Number of raw single end reads:30465	
#Number of mapped paired reads:77001	
#Number of mapped unpaired reads:11003	
#Number of mapped single end reads:	30465
#Number of mapped paired reads pcr duplicates:	18922
#Number of mapped unpaired reads pcr duplicates:	1726
#Number of mapped single end reads pcr duplicates:	2471
#Number of ends at mapping quality 0-255:	188413:0	36:1	 ...

vals = []
dist_data = {}
$stdin.readlines.each do |l| 
  if l =~ /reads/
    k, v = l.split(":")
    vals << v.gsub(/\s+/, '')
  end
  if l =~ /mapping quality/
    p = l.split
    p.pop; p.pop
    p.each do |e|
      if e =~ /^\d/ && e != "0-255:"
        amount, mapping_qual = e.split(":")
        dist_data[mapping_qual]=amount
      end
    end
  end
end

json_map_stats = <<EOF
{
  "raw_reads"             : #{vals[0]},
  "raw_single"            : #{vals[1]},
  "mapped_paired"         : #{vals[2]},
  "paired_one_end_mapped" : #{vals[3]},
  "mapped_single"         : #{vals[4]},
  "paired_dups"           : #{vals[5]},
  "unpaired_dups"         : #{vals[6]},
  "single_dups"           : #{vals[7]}
}
EOF

json_data = []
dist_data.each do |mapq, amount|
  json_data << " { \"x\":#{mapq}, \"y\":#{amount} }"
end
json_dist_mapq = "[\n" + json_data.join(",") + "]\n"

File.open("map_stats.json", "w") {|f| f.write(json_map_stats)}
File.open("dist_mapq.json", "w") {|f| f.write(json_dist_mapq)}
