#!/bin/bash
#
set -e

#bw_path="/tmp/drd/breakway.0.7"
bw_path="/data/rogers/drio_scratch/bb/local/breakway.0.7"
dt_file="/tmp/drd/dtrans.txt"
bam="/data/rogers2/drio_scratch/gibbon/AWG_NLEU.00_000pA/AWG_NLEU.00_000pA.merged.bam"
ff="/data/rogers2/drio_scratch/genomes/Nleu1.0.fixed.fa"
fai="/data/rogers2/drio_scratch/genomes/Nleu1.0.fixed.fa.fai"
expect_window=4000
dt_inc=100
dt_bin=500
rl=50 # read length
g_size="2936035333" # genome size in bp

log()
{
  echo "`date` >> $1" >&2
}

testing_reads()
{
  local few_reads="100000"
  log "WARNING: testing mode"
  log "extracting few alignments ($few_reads)"
  samtools view -h $bam | head -${few_reads} | samtools view -bS - > i.bam 2>/dev/null
  bam="i.bam"
}

find_ped_stats()
{
  log "finding PED stats [./pen.txt]"
  perl $bw_path/scripts/breakway.parameters.pl 0 $expect_window $bam > pen.txt 2> /dev/null
  ped_mean=`cat ped.txt  | grep Mean     | awk  -F= '{print $2}' | sed 's/ //g'`
  ped_std=`cat ped.txt   | grep Standard | awk  -F= '{print $2}' | sed 's/ //g'`
  ped_lower=`cat ped.txt | grep lower    | awk  -F= '{print $2}' | sed 's/ //g'`
  ped_upper=`cat ped.txt | grep upper    | awk  -F= '{print $2}' | sed 's/ //g'`
  ped_lower=$[$lower+100]
  ped_upper=$[$upper+100]
  log "PED mean: ${ped_mean}" 
  log "PED std : ${ped_std}"
}

get_bam_stats()
{
  log "bam stats [./stats.txt]"
  dbamstats $bam > stats.txt 2>/dev/null
  nm_pe_reads=`cat stats.txt | grep "mapped paired reads:" | awk -F: '{print $2}'`
  min_cs=`echo "$ped_lower/$dt_inc" | bc`
  max_cs=`echo "$ped_upper/$dt_inc" | bc`
  hap_ccov=`echo "scale=2;(($nm_pe_reads*(2*$rl+$ped_mean)) / $g_size) / 2" | bc`
  std_ccov=`echo "$nm_pe_reads*$ped_std/2" | bc`
  log "# mapped pe reads      : ${nm_pe_reads}"
  log "haploid coverage       : ${hap_ccov}"
  log "st dev haploid coverage: ${std_ccov}"
}

run_dtranslocations()
{
  log "running dtranslocations [./dtrans.txt]"
  dtranslocations -i $dt_inc -b $dt_bin -l $ped_lower -L $ped_upper $bam > dtrans.txt 2>/dev/null
}

find_genome_size()
{
  log "finding genome size"
  #g_size=`cat $ff | ruby -ne 'BEGIN{@t=0}; @t = @t + $_.chomp.size unless $_ =~ /^>/; END{puts @t.to_i}'` 
  log "$g_size bp"
}

find_svs()
{
  log "finding sv [./breakway.out.bw]"
  #cat <<-EOF
  perl $bw_path/breakway.run.pl \
  --dtransfile $dt_file \
  --mincs $min_cs --maxcs $max_cs \
  --interval $dt_inc \
  --cutoff 4 \
  --bamfile $bam \
  --fai $fai \
  --filehandle breakway \
  --bwfolder $bw_path \
  --mindist $ped_lower --maxdist $ped_upper --mean $hap_ccov --stdev $std_ccov --score 0.01 \
  --strand SAME --ordering 21 \
  > breakway.out.bw 2> /dev/null
  #EOF
  #--segdupsfile [genomicSuperDups.txt] --selfchainfile [chainSelf.txt] --repmaskerfile [rmsk.txt] 
}

#
# Main
testing_reads
find_ped_stats
get_bam_stats
run_dtranslocations
find_genome_size
find_svs
log "done"
