#!/usr/bin/env perl

use strict;
use warnings;


my $usage = qq{
Usage: bfast2novo.pl <in.fastq> <out.prefix> <split num>

Note: <in.fastq> must have been generated from the solid2fastq 
utility in the BFAST distribution.  <split num> can be omitted,
resulting in one large file.

};

die($usage) if (@ARGV != 2 && @ARGV != 3);
my ($in_fastq, $out_prefix, $split_num);
my $warn = 0;
if(@ARGV == 2) {
	($in_fastq, $out_prefix) = @ARGV;
	$split_num = 0;
}
else { # must equal 3
	($in_fastq, $out_prefix, $split_num) = @ARGV;
	if($split_num <= 0) {
		die("<split num> was <= 0 [$split_num]\n");
	}
}

my ($pe_n, $se_n) = (0, 0);
my ($pe_i, $se_i) = (1, 1);
my $ctr=0;

open(FHin, "$in_fastq") || die("** Could not open '$in_fastq'.\n");

if(0 == $split_num) {
	open(FHout_pe1, ">$out_prefix.read1.fastq") || die;
	open(FHout_pe2, ">$out_prefix.read2.fastq") || die;
	open(FHout_se, ">$out_prefix.single.fastq") || die;
}
else {
	open(FHout_pe1, ">$out_prefix.read1.$pe_i.fastq") || die;
	open(FHout_pe2, ">$out_prefix.read2.$pe_i.fastq") || die;
	open(FHout_se, ">$out_prefix.single.$se_i.fastq") || die;
}

my (@read1, @read2);
@read1 = &read_fastq(\*FHin);
@read2 = &read_fastq(\*FHin);

while(@read1 && @read2) {
	if($read1[0] eq $read2[0]) { # two ends
		if(0 != $split_num && $split_num <= $pe_n) {
			$pe_i++;
			$pe_n=0;
			close(FHout_pe1); close(FHout_pe2);
			open(FHout_pe1, ">$out_prefix.read1.$pe_i.fastq") || die;
			open(FHout_pe2, ">$out_prefix.read2.$pe_i.fastq") || die;
		}
		print FHout_pe1 $read1[1];
		print FHout_pe2 $read2[1];
		$pe_n++;
		@read1 = &read_fastq(\*FHin);
		@read2 = &read_fastq(\*FHin);
	}
	else {
		if(0 != $split_num && $split_num <= $se_n) {
			$se_i++;
			$se_n=0;
			open(FHout_se, ">$out_prefix.single.$se_i.fastq") || die;
		}
		if($read1[0] le $read2[0]) {
			print FHout_se $read1[1];
			$se_n++;
			@read1 = &read_fastq(\*FHin);
		} 
		else {
			print FHout_se $read2[1];
			$se_n++;
			@read2 = &read_fastq(\*FHin);
		}
	}
}
while(@read1) {
	if(0 != $split_num && $split_num <= $se_n) {
		$se_i++;
		$se_n=0;
		open(FHout_se, ">$out_prefix.single.$se_i.fastq") || die;
	}
	print FHout_se $read1[1];
	$se_n++;
	@read1 = &read_fastq(\*FHin);
}
while(@read2) {
	if(0 != $split_num && $split_num <= $se_n) {
		$se_i++;
		$se_n=0;
		open(FHout_se, ">$out_prefix.single.$se_i.fastq") || die;
	}
	print FHout_se $read2[1];
	$se_n++;
	@read2 = &read_fastq(\*FHin);
}
close(FHout_pe1);
close(FHout_pe2);
close(FHout_se);

sub read_fastq {
	my $fh = shift;
	my ($name, $seq, $comment, $qual);
	if(!defined($name = <$fh>) 
		|| !defined($seq = <$fh>) 
		|| !defined($comment = <$fh>) 
		|| !defined($qual = <$fh>)) {
		return ();
	}
	chomp($name); chomp($seq); chomp($comment); chomp($qual);
	my $read = qq/$name\n$seq\n$comment\n$qual\n/;
	#if($name =~ m/^@(\d+)_(\d+)_(\d+)/) {
	if($name =~ m/^@\D+(\d+)_(\d+)_(\d+)/) {
		my $key = sprintf("%.4d_%.4d_%.4d", $1, $2, $3); # this line could be improved on 64-bit machines
		return ($key, $read);
	}
	else {
		my $key = 0; # not recognized
		if(0 == $warn) {
			print STDERR "Warning: read name could not be parsed:\n";
			print STDERR "Name [$name]\n";
			print STDERR "Subsequent warnings will be surpressed.\n";
			$warn = 1;
		}
		return ($key, $read);
	}
}
