#!/usr/bin/env ruby19
#
# time samtools view ./0.bam | ./plot_dist_mapq.rb
# for i in `seq 0 8`; do samtools view ./$i.bam | ./plot_dist_mapq.rb $i  2>/dev/null & done ; wait &
#
# /stornext/snfs1/next-gen/software/jdk1.6.0_01/bin/java -jar -Xmx8g 
# /stornext/snfs1/next-gen/software/picard-tools/current/MergeSamFiles.jar 
# INPUT=1.bam INPUT=2.bam
# TMP_DIR=/space1/tmp/ OUTPUT=./out.bam
# VALIDATION_STRINGENCY=STRICT
#
def per(x,t)
  p = (x * 100)/t
  "#{p}%"
end

def usage
  puts <<END
Usage: samtools view ./foo.bam | #{__FILE__} <seed_output_file>
END
  exit 1
end

ARGV.size != 1 && usage
output_seed = ARGV[0]

data_fn  = "data_#{output_seed}.txt"
cfg_fn   = "cfg_#{output_seed}.txt" 
gout_fn  = "#{output_seed}.png" 
mq_fn    = "map_qual_#{output_seed}.txt" 
u        = 0

h = Hash.new(0)
t = 0
STDIN.each do |l|
  $stderr.printf("\r%d", t) if t % 50000 == 0
  flag  = l.split[1].to_i
  unless flag == 4
    mq    = l.split[4].to_i
    h[mq] = h[mq] + 1
  else
    u = u + 1
  end 
  t = t + 1
end

m = t - u
mq_stats =<<EOF
-----------------------
Total   : #{t}
mapped  : #{m}\t #{per(m,t)}
Unmapped: #{u}\t #{per(u,t)}
0  MAPQ : #{h[0]}\t #{per(h[0], m)}
!0 MAPQ : #{m - h[0]}\t #{per(m-h[0], m)}
-----------------------
EOF
File.open(mq_fn, "w") {|f| f.puts mq_stats }

File.open(data_fn, "w") do |df|
  h.sort{|a,b| a[0]<=>b[0]}.each do |v|
    df.puts "#{v[0]}\t#{v[1]}"
  end
end

File.open(cfg_fn, "w") do |cfg| 
  cfg.puts DATA.read.gsub(/XXOUTXX/, gout_fn).gsub(/XXDATAXX/, data_fn)
end

`cat ./#{cfg_fn} | gnuplot`
`rm -f #{cfg_fn}`

__END__
set terminal png nocrop size 800,300
#set key left top
set style fill solid border -1
set boxwidth 1.0
set nokey
set xtics 10 
#set grid x y
#set datafile separator ","
set t png
set log y
set xrange [-1:257]
#set yrange [0:100]
set output "XXOUTXX"
set ylabel "# reads (alignments) -- log scale"
set xlabel "Mapping quality"
set title "Mapping quality dist. (XXOUTXX)"
plot "XXDATAXX" notitle with boxes
