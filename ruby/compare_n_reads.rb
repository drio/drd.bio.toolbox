#!/usr/bin/env ruby
#
require 'json'

def usage(msg=nil)
  $stderr.puts "ERROR: #{msg}" if msg
  $stderr.puts "Usage: cat recipe.json | #{$0} merged_bam"
  exit 1 if msg
end

def l(msg)
  $stderr.puts "#{Time.now}>> #{msg}"
end

def find_n_reads(b)
  cmd = "samtools view #{b} | wc -l" 
  l "Iterating over #{b}"
  begin
    n_reads = `#{cmd}`
  rescue
    $stderr.puts "Problems using samtools to extract n_reads"; raise
  end
  n_reads.to_i
end

usage("Wrong # of params") unless ARGV.size == 1
merged_bam = ARGV[0]
usage("I can't open file: #{merged_bam} ") unless File.exists? merged_bam
begin
  j = JSON.parse($stdin.read)
rescue
  $stderr.puts "Problems parsing JSON from stdin"; raise
end

count = 0
j['bams'].each do |b| 
  count += find_n_reads b
  l "count = #{count}"
end
m_count = find_n_reads merged_bam
l "count = #{m_count}"

puts "#_reads_single_bams, #_reads_merged_bam, 1-2"
puts "#{count}, #{m_count}, #{count - m_count}"

__END__
{
  "ref_fasta"   : "/stornext/snfs0/rogers/drio_scratch/genomes/marmoset.bwa.pipe.fixed.fa",
  "bams"        : [
    "/stornext/snfs0/rogers/drio_scratch/marmoset/illumina/samples/32789/IWG_MDPJR.32789-1_1pA.l7.4.fixed.bam",
    "/stornext/snfs0/rogers/drio_scratch/marmoset/illumina/samples/32789/IWG_MDPJR.32789-1_1pA.l5.2.fixed.bam",
    "/stornext/snfs0/rogers/drio_scratch/marmoset/illumina/samples/32789/IWG_MDPJR.32789-1_1pA.l4.1.fixed.bam",
    "/stornext/snfs0/rogers/drio_scratch/marmoset/illumina/samples/32789/IWG_MDPJR.32789-1_1pA.l6.3.fixed.bam"],
  "title"       : "bn.marmoset.32789",
  "prj_name"    : "marmoset",
  "sample_name" : "32789"
}

