#!/usr/bin/env ruby

require File.expand_path(File.dirname(__FILE__)) + "/common"

include Help
include Common
usage=<<-EOF
Usage: #{File.basename $0} <ref_path> <prj_name> <sample_name>
EOF
Help.set_usage_text usage

Help.error "Wrong number of params" unless ARGV.size == 3 
ref, prj, sample = ARGV
Help.error "Ref path not found" unless File.exists?(ref)
Help.error "No data in stdin." unless data_in_stdin?

class Array
  def to_json
    '[' + "\n" + '"' + self.join('",' + "\n" +'"') +  '"]'
  end
end

bams=[]
$stdin.each_line do |b|
  bams << File.expand_path(b.chomp)
end

puts DATA.read
     .gsub(/REF/   , ref)
     .gsub(/PRJ/   , prj)
     .gsub(/SAMPLE/, sample)
     .gsub(/BAMS/  , bams.to_json)

__END__
{
  "ref_fasta"   : "REF",
  "bams"        : BAMS,
  "title"       : "bn.PRJ.SAMPLE",
  "prj_name"    : "PRJ",
  "sample_name" : "SAMPLE"}
}
