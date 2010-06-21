#!/usr/bin/env ruby19
#
# run_bwa.rb: wrapper on top of bwa
#
# /stornext/snfs1/next-gen/software/bwa-svn/solid2fastq.pl 0044_20100421_2_SP_ANG_AUT_N_200716041_1_1sA_01003310991_3_ bwa
# gzip -d bwa.read1.fastq.gz &
# gzip -d bwa.read2.fastq.gz &
# wait
# mkdir -p splits.read1
# mkdir -p splits.read2
# split -a 1 -d -l 40000000 ../bwa.read1.fastq ./splits.read1/read1.fastq. &
# split -a 1 -d -l 40000000 ../bwa.read2.fastq ./splits.read2/read2.fastq. &
# wait
# bwa aln -c -t8 $ref $fq > $i.sai
# bwa samse $ref $i.sai $fq > $i.sam
# samtools view -bS $i.sam > $i.bam
# samtools flagstat $i.bam  > $i.stats.txt
#
require 'optparse'
require 'ostruct'

class String
  def add_n(s, nn=1)
    ns = "\n" * nn
    self.replace("#{self}#{s}#{ns}")
  end
end

class BWA_wrapper
  def initialize(arguments)
    @options   = OpenStruct.new
    @arguments = arguments
    @s = "".add_n "#!/bin/bash"
    @s.add_n "#"
    @s.add_n "set -e", 2
    check_binaries
    set_toolbox
    log "Parsing options... "
    exit 1 unless parsed_options?
    print_values
  end

  def run
    log "Arguments look good."
    generate_script
    puts @s
  end

  def parsed_options?
    # Specify options
    opts = OptionParser.new
    opts.on('-h', '--help')         { output_help }
    opts.on('-r', '--path_raw r')   {|r| @options.path_raw = r }
    opts.on('-u', '--r_name   u')   {|u| @options.run_name = u }
    opts.on('-s', '--s2f      s')   {|s| @options.s2f = s }
    opts.on('-b', '--bb       b')   {|b| @options.bb  = b }
    opts.on('-f', '--ref      f')   {|f| @options.ref = f }
    
    opts.parse!(@arguments) rescue return false 
    process_options
    true  
  end   

  private

  def check_binaries
    %w(samtools bsub).each do |b|
      p = `which #{b}`
      log "Can't find bin: #{b}", 1, 1 unless $? == 0
      log "#{b} path: #{p}"
    end
  end

  def set_toolbox
    @tb_dir = File.dirname(__FILE__).gsub(/\//,'')
    log "Toolbox path: #{@tb_dir}"
    %w(bam_sort_dups.sh bam.stats.sh plot_dist_mapq.rb).each do |f|
      log "Toolbox, item not found {#{f}}", 1, 1 unless File.executable?(@tb_dir + "/" + f)
    end
  end

  def generate_script
    fastqs
    ungzip
    aln
    sam
    bam_sort_dups
    stats
    mq_plot
  end

  def fastqs
    @s.add_n "#{@s2f} #{@run_name} bwa"
  end

  def ungzip
    @s.add_n ""
    @s.add_n "gzip -d bwa.read1.fastq.gz &"
    @s.add_n "gzip -d bwa.read2.fastq.gz &"
    @s.add_n "wait"
  end

  def aln
    @s.add_n ""
    @s.add_n "#{@bb} aln -c -t8 #{@ref} bwa.read1.fastq > read1.sai"
    @s.add_n "#{@bb} aln -c -t8 #{@ref} bwa.read2.fastq > read2.sai"
  end

  def sam
    @s.add_n ""
    @s.add_n "#{@bb} sampe #{@ref} read1.sai read2.sai bwa.read1.fastq bwa.read2.fastq > #{fixed_run_name}.sam"
  end

  # raw/0097_20100111_1_SP_ANG_RC_LCA_K86_3_1sA_01003280864_2_.sam
  def fixed_run_name
    (File.basename @run_name).gsub(/_$/, '')
  end

  def bam_sort_dups
    @s.add_n ""
    @s.add_n "#{@tb_dir}/bam_sort_dups.sh #{fixed_run_name}.sort.dups #{fixed_run_name}.sam | /bin/bash"
  end

  def stats
    @s.add_n ""
    @s.add_n "#{@tb_dir}/bam.stats.sh -i #{fixed_run_name}.sort.dups.bam -e 2 -m 1 > #{fixed_run_name}.stats.txt"
  end

  def mq_plot
    b = fixed_run_name + ".sort.dups.bam"
    @s.add_n ""
    @s.add_n "samtools view #{b }| #{@tb_dir}/plot_dist_mapq.rb #{b}"
  end

  def process_options
    @path_raw = @options.path_raw || nil
    @run_name = @options.run_name || nil
    @s2f      = @options.s2f      || nil
    @bb       = @options.bb       || nil
    @ref      = @options.ref      || nil
    log "Problems with path to raw data"    , 1, 1 if @path_raw.nil? || !Dir.exists?(@path_raw)
    log "No run name provided"              , 1, 1 if @run_name.nil?
    log "Problems using solid2fastq script" , 1, 1 if @s2f.nil?      || !File.executable?(@s2f)
    log "Problems with bwa bin"             , 1, 1 if @bb.nil?       || !File.executable?(@bb)
    log "Problems with ref genome"          , 1, 1 if @ref.nil?      || !File.exists?(@ref.gsub(/\.cs/,'') + ".fa")
  end

  def output_help
    puts DATA.read
  end

  def print_values
    log "raw path        : #{@path_raw}"
    log "run name        : #{@run_name}"
    log "s2f             : #{@s2f}"
    log "bwa bin         : #{@bb}"
    log "ref genome seed : #{@ref}"
  end
  
  def log(msg, bye=0, help=0)
    $stderr.puts "LOG: " + msg 
    output_help   if help != 0
    exit bye.to_i if bye  != 0
  end
end

BWA_wrapper.new(ARGV).run

__END__

Usage:
  run_bwa.rb [options]

Options:
 -h, --help          Displays help message

 -r, --path_raw r    Path to raw data
 -u, --run_name u    Run name
 -s, --s2f s         Path to solid2fastq script 
 -b, --bb  b         Path to bwa bin
 -r, --ref r         Path to ref genome

Examples:

$ # Interactive
$ run_bwa.rb -r ./raw -u "raw/0097_20100111_1_SP_ANG_RC_LCA_K86_3_1sA_01003280864_2_" -s /stornext/snfs1/next-gen/software/bwa-svn/solid2fastq.pl -b /stornext/snfs1/next-gen/software/bwa-svn/bwa -f /stornext/snfs4/next-gen/solid/bwa.references/h/hsap.36.1.hg18/hsap_36.1_hg18.cs | /bin/bash"

$ # lsf mode
$ bsub -q high -J great_times "../run_bwa.rb -r ./raw -u "raw/0097_20100111_1_SP_ANG_RC_LCA_K86_3_1sA_01003280864_2_" -s /stornext/snfs1/next-gen/software/bwa-svn/solid2fastq.pl -b /stornext/snfs1/next-gen/software/bwa-svn/bwa -f /stornext/snfs4/next-gen/solid/bwa.references/h/hsap.36.1.hg18/hsap_36.1_hg18.cs | /bin/bash"
