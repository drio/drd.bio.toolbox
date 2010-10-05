#!/bin/env ruby

require File.expand_path(File.dirname(__FILE__)) + "/common"

def usage
  puts "Usage: #{__FILE__} eval1.txt ... eval_n.txt"
  exit 1
end

hash_common = Common::load_common
ARGV.size < 1 && usage

total = 0
accuracy_data = {}
ARGV.each {|fn| accuracy_data[fn] = {} }
ARGV.each do |fn| 
  File.open(fn).each_line do |l|
    l_data = l.split
    accuracy_data[fn][:mapq]   = l_data[0]
    accuracy_data[fn][:ok]     = l_data[4]
    accuracy_data[fn][:mapped] = l_data[5]
    total = l_data[6] if total == 0
  end
end

# Header
printf "mapq "
ARGV.each {|fn| printf "OK_#{fn} mapped_#{fn} " }
printf "total\n"

# Data
(0..255).inject([]) {|set, n| set << n }.reverse.each do |mapq|
  printf "#{mapq} "
  ARGV.each {|fn| printf "#{accuracy_data[fn][:ok]} #{accuracy_data[fn][:mapped]} "}
  printf "#{total}\n"
end
