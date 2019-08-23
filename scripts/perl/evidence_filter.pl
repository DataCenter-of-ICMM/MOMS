#!/usr/bin/env perl
use strict;
use warnings;
use POSIX;

if(@ARGV < 3){
	print "Usage: $0 evidences.txt cnt ratio [maxdepth]\n";
	exit 1;
}

my $evidencefile = $ARGV[0];
my $count_threshold = $ARGV[1];
my $ratio_threshold = $ARGV[2];
if($count_threshold !~ /^[0-9]+$/){
	print STDERR "ERROR: the count threshold must be a number.\n";
	exit 1;
}
if($count_threshold < 3){
	$count_threshold = 3;
}
my $count_threshold2 = $count_threshold * 2;
if($ratio_threshold !~ /^[0-9.]+$/ or $ratio_threshold > 1){
	print STDERR "ERROR: the ratio threshold must be a number between 0 and 1.\n";
	exit 1;
}
my $max_depth = 9 ** 99;
if(@ARGV >= 4){
	if($ARGV[3] !~ /^[0-9.]+$/){
		print STDERR "ERROR: the maximum depth must be a number.\n";
		exit 1;
	}
	$max_depth = $ARGV[3];
	if($max_depth < $count_threshold2 * 2){
		print STDERR "ERROR: the maximum depth must be no less than 4 times of the count threshold.\n";
		exit 1;
	}
}

my ($key, $st, $new_indices);
my ($last_ctg1, $last_ctg2, $last_end1, $last_end2);
my ($ctg1, $ctg2, $end1, $end2);
my ($line, @columns);
my @hits = ();
my @local_hits = ();
my @distances = ();
my %stats = ();
# to filter outliers
open(IN, "<$evidencefile") or die("Can not open $evidencefile for reading\n");
while($line = <IN>){
	chomp($line);
	if($line =~ /^#/){
		print "$line\n";
		next;
	}
	@columns = split(/\t/, $line);
	next if(scalar(@columns) < 9);
	($ctg1, $ctg2, $end1, $end2) = @columns[1..4];
	if(defined $last_ctg1 && ($ctg1 != $last_ctg1 || $ctg2 != $last_ctg2 || $end1 ne $last_end1 || $end2 ne $last_end2)){
		$key = "${last_ctg1}_${last_ctg2}_${last_end1}_${last_end2}";
		$st = calculateStatistics(\@distances);
		$new_indices = filter(\@distances, $st);
		@local_hits = @local_hits[@$new_indices];
		@distances = @distances[@$new_indices];
		$st = calculateStatistics(\@distances);
		if( ($st->{mean} > -3000) && ($st->{std} < 10000) ){
			push(@hits, @local_hits);
		}
		@local_hits = ();
		@distances = ();
	}
	($last_ctg1, $last_ctg2) = ($ctg1, $ctg2);
	($last_end1, $last_end2) = ($end1, $end2);
	push(@distances, $columns[8]); # distance
	push(@local_hits, $line);
}
close IN;
if(defined $last_ctg1){
	$key = "${last_ctg1}_${last_ctg2}_${last_end1}_${last_end2}";
	$st = calculateStatistics(\@distances);
	$new_indices = filter(\@distances, $st);
	@local_hits = @local_hits[@$new_indices];
	@distances = @distances[@$new_indices];
	$st = calculateStatistics(\@distances);
	if( ($st->{mean} > -3000) && ($st->{std} < 10000) ){
		push(@hits, @local_hits);
	}
}

# to filter evidences based on counts and consistency
my %indices = ("5'5'" => 0, "5'3'" => 1, "3'5'" => 2, "3'3'" => 3);
my (@counts, $total_count);
@counts = ((0) x 4);
$total_count = 0;
($last_ctg1, $last_ctg2) = (undef, undef);
foreach $line (@hits){
	@columns = split(/\t/, $line);
	next if(scalar(@columns) < 9);
	($ctg1, $ctg2, $end1, $end2) = @columns[1..4];
	if(defined $last_ctg1 && ($ctg1 != $last_ctg1 || $ctg2 != $last_ctg2) ){
		if($total_count >= $count_threshold && $total_count <= $max_depth){
			@counts = sort {$b <=> $a} (@counts);
			if( $counts[0] == $total_count ||
				($total_count >= $count_threshold2) && ($counts[0] >= $total_count * $ratio_threshold) ){
				foreach (@local_hits){
					print "$_\n";
				}
			}
		}
		@counts = ((0) x 4);
		$total_count = 0;
		@local_hits = ();
	}
	($last_ctg1, $last_ctg2) = ($ctg1, $ctg2);
	my $idx = $indices{"$end1$end2"};
	$counts[$idx]++;
	$total_count++;
	push(@local_hits, $line);
}
if(defined $last_ctg1){
	if($total_count >= $count_threshold && $total_count <= $max_depth){
		@counts = sort {$b <=> $a} (@counts);
		if( $counts[0] == $total_count ||
			($total_count >= $count_threshold2) && ($counts[0] >= $total_count * $ratio_threshold) ){
			foreach (@local_hits){
				print "$_\n";
			}
		}
	}
}
exit 0;

sub filter
{
	my ($distances, $stats) = @_;
	if(scalar(@$distances) <= 3){
		return [0..$#distances];
	}
	my @new_distances = ();
	my ($mean, $std, $cnt) = ($stats->{mean}, $stats->{std}, $stats->{cnt});
	my @diffs = ();
	foreach (@$distances){
		push(@diffs, abs($_ - $mean));
	}
	my @sorted = sort{$diffs[$a] <=> $diffs[$b]} (0..$#distances);
	my $is = $#sorted - ceil(scalar(@sorted) * 0.1);
	my $nExtended = 0;
	while($is < scalar(@sorted)){
		if( $diffs[$sorted[$is]] >= 1.2 * $std ){
			if($nExtended > 0){
				$is--;
			}
			last;
		}
		last if($is == $#sorted);
		$is++;
		$nExtended++;
	}
	return [sort @sorted[0..$is]];
}

sub calculateStatistics
{
	my ($distances) = @_;
	my @sorted = ();
	my ($mean, $std);
	my $sum = 0;
	foreach (@$distances){
		$sum += $_;
	}
	$mean = sprintf("%.1f", $sum / scalar(@$distances));
	$sum = 0;
	foreach (@$distances){
		$sum += ($_ - $mean) * ($_ - $mean);
	}
	$std = sprintf("%.1f", sqrt($sum / scalar(@$distances)));

	return {mean=>$mean, std=>$std, cnt=>scalar(@$distances)};
}

