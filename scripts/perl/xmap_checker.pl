#!/usr/bin/env perl
# Copyright (C) 2018,2019 Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences
# MOMS is licensed under the Mulan PSL v1.
# You can use this software according to the terms and conditions of the Mulan PSL v1.
# You may obtain a copy of Mulan PSL v1 at:
#    http://license.coscl.org.cn/MulanPSL
# THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT, MERCHANTABILITY OR FIT FOR A PARTICULAR
# PURPOSE.
# See the Mulan PSL v1 for more details.
use strict;

my $usage = << "USAGE";
Usage: $0 sample.xmap
USAGE

if(scalar(@ARGV) ne 1){
	print $usage;
	exit 1;
}

my $filename = shift;

open(XM,"<$filename") or die("Can not open \"$filename\" for reading\n");

my $line;
my @columns = ();
my %qry = ();
my %ref = ();
while($line = <XM>){
	chomp($line);
	next if($line =~ /^#/);
	@columns = split(/\t/, $line);
	next if(scalar(@columns) < 3);
	$qry{$columns[1]} = 1;
	$ref{$columns[2]} = 1;
}
close XM;

printf "%5s BNG contigs have been flagged as conflicting\n", scalar(keys %qry);
printf "%5s NGS contigs have been flagged as conflicting\n", scalar(keys %ref);

exit 0;
