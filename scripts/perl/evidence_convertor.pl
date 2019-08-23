#!/usr/bin/env perl
use strict;
use warnings;

if(@ARGV < 3){
	print "Usage: $0 keys.txt evidences.txt insert_size\n";
	exit 1;
}

my $keyfile = $ARGV[0];
my $evidencefile = $ARGV[1];
my $insert_size = $ARGV[2];

my %lengths = ();
my $line;
my @columns = ();
my ($id, $name, $len);
open(IN, "<$keyfile") or die("Can not open $keyfile for reading\n");
while($line = <IN>){
	chomp($line);
	next if($line !~ /^\d/);
	@columns = split(/\t/, $line);
	next if(scalar(@columns) < 3);
	($id, $name, $len) = @columns[0..2];
	$lengths{$id} = $len;
}
close IN;

my ($ctg1, $ctg2, $end1, $end2, $dist);
my ($len1, $len2, $adj1, $adj2, $span);
my ($s1, $e1, $s2, $e2);
open(IN, "<$evidencefile") or die("Can not open $evidencefile for reading\n");
my $rlen = 100;
while($line = <IN>){
	chomp($line);
	next if($line =~ /^#/);
	@columns = split(/\t/, $line);
	next if(scalar(@columns) < 8);
	(undef, $ctg1, $ctg2, $end1, $end2, undef, undef, $dist) = @columns;
	($len1, $len2) = ($lengths{$ctg1}, $lengths{$ctg2});
	$span = int($insert_size - $dist);
	$adj1 = sprintf("%d", int($len1/($len1+$len2)*$span));
	$adj2 = $span - $adj1;
	$len1 -= $adj1;
	$len2 -= $adj2;
	if($end1 eq "5'"){
		if($end2 eq "5'"){ # 5'5'
			($s1, $e1, $s2, $e2) = ($len1, $len1-$rlen, $len2, $len2-$rlen);
		}
		else{ # 5'3'
			($s1, $e1, $s2, $e2) = ($len1, $len1-$rlen, $len2-$rlen, $len2);
		}
	}
	else{
		if($end2 eq "5'"){ # 3'5'
			($s1, $e1, $s2, $e2) = ($len1-$rlen, $len1, $len2, $len2-$rlen);
		}
		else{ # 3'3'
			($s1, $e1, $s2, $e2) = ($len1-$rlen, $len1, $len2-$rlen, $len2);
		}
	}
	print "$ctg1\t$s1\t$e1\t$ctg2\t$s2\t$e2\n";
}
close IN;

exit 0;

