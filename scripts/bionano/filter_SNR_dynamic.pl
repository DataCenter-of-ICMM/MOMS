#!/usr/bin/perl -w
#############################################################################
# Last modified date: 08/08/2017                                            #
# Purpose: Molecule filtering using dynamically estimated SNR for BNX files #
#                                                                           #
# Author: Xiang Zhou, Computational Biologist                               #
# Email : xzhou@bionanogenomics.com                                         #
# Affiliation: Research Department, BioNano Genomics Inc.                   #
#                                                                           #
# Usage:                                                                    #
#   perl filter_SNR_dynamic.pl [options] <Args>                             #
# Options:                                                                  #
#   -h    : This help message                                               #
#   -i    : Input BNX file (Required, can be 1 or 2 colors)                 #
#   -o    : Output BNX file (default: does not generate new BNX)            #
#   -C    : Only filter the BNX on the channel specified (use with -o)      #
#   -K    : Output BNX file BUT KEEP ALL labels (default: OFF, valid w/ -o) #
#   -m    : Minimum allowed dynamic SNR cutoff (default: 2)                 #
#   -M    : Maximum allowed dynamic SNR cutoff (default: 20)                #
#   -d    : Default dynamic SNR cutoff (3.5) if outside the range of [m, M] #
#   -L    : Minimum number of labels for the script to estimate SNR         #
#   -P    : Plot SNR histogram using the PDF filename specified (no default)#
#############################################################################

use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Spec;
use File::Path qw(make_path);
use Cwd qw(abs_path realpath);
use Getopt::Std;

my $progpath = abs_path(dirname($0));

my ($default_min_SNR, $default_max_SNR, $default_SNR, $min_num_of_labels, $filter_channel) = (2, 20, 3.5, 0, 0);
my ($header, $command, $RefAligner) = ("") x 3;
my ($fin, $fout, $fout_summary, $fout_SNR);
my %opts;

Init();
Check_inputs();

if( defined($opts{o}) && ! $opts{K} ){
	$RefAligner = "$progpath/binary/RefAligner";
	if(! -f $RefAligner){
		print STDERR ("ERROR: Could not find RefAligner in the default path!\n");
		exit (1);
	}
}

######################################################
my (@OUT, @SNR_cutoff);
my ($IN, $OUT, $OUT_summary);
my ($channel_1, $channel_2);
my ($count, $total, $remain) = (1, 0, 0);
my (%mol, %label, %SNR, %mol_SNR);
######################################################
# Dynamic SNR cutoff estimation
if( defined($fout_SNR) ){
	@SNR_cutoff = GetSNR($fin, $fout_SNR);
}
else{
	@SNR_cutoff = GetSNR($fin);
}

# Static SNR cutoff
# @SNR_cutoff = (YOUR_SNR_CUTOFF_VALUEs);

#print "Dynamic SNR cutoff: ", join(" ", @SNR_cutoff), "\n";
if( $SNR_cutoff[0] eq "NA" || $SNR_cutoff[0] != -1 && ($SNR_cutoff[0] < $default_min_SNR || $SNR_cutoff[0] > $default_max_SNR) ){
	print STDERR ("Warning: The estimated SNR cutoff value for channel 1 is out of range! Use default SNR cutoff value of $default_SNR instead.\n");
	$SNR_cutoff[0] = $default_SNR;
}
if ($SNR_cutoff[1] eq "NA") {
    splice @SNR_cutoff, 1;
}
else{
	if( $SNR_cutoff[1] != -1 && ($SNR_cutoff[1] < $default_min_SNR || $SNR_cutoff[1] > $default_max_SNR) ){
		print STDERR ("Warning: The estimated SNR cutoff value for channel 2 is out of range! Use default SNR cutoff value of $default_SNR instead.\n");
		$SNR_cutoff[1] = $default_SNR;
	}
}
######################################################
#open($OUT_summary, ">".$fout_summary) || die ("ERROR: Can't open $fout_summary: $!\n");
print $OUT_summary("PROGRAM PARAMETERS\n");
print $OUT_summary "Input BNX file:\t", $fin, "\n";
if(defined($opts{o})){
	print $OUT_summary "Output BNX file:\t", $fout, "\n";
}
print $OUT_summary "Minimum allowed dynamic SNR cutoff:\t", $default_min_SNR, "\n";
print $OUT_summary "Maximum allowed dynamic SNR cutoff:\t", $default_max_SNR, "\n";

for(my $i = 1; $i <= @SNR_cutoff; $i++){
	if($SNR_cutoff[$i-1] == -1){
		print $OUT_summary "Estimated SNR cutoff for channel $i:\t", -1, "\n";
	}
	else{
		print $OUT_summary "Estimated SNR cutoff for channel $i:\t", sprintf("%.3f", $SNR_cutoff[$i-1]), "\n";
	}
}
if(! defined($opts{o})){
	exit 0;
}
my $tmp = sprintf("%.3f", $SNR_cutoff[0]);
if($opts{K}){
	for(my $i = 1; $i <= @SNR_cutoff; $i++){
		# print $OUT_summary "Estimated SNR cutoff for channel $i:\t", sprintf("%.3f", $SNR_cutoff[$i-1]), "\n";
		my $tmp;
		if($SNR_cutoff[$i-1] == -1){
			$tmp = -1;
		}
		else{
			$tmp = sprintf("%.3f", $SNR_cutoff[$i-1]);
		}
		$command = "echo \"# Estimated SNR cutoff value for channel $i using histogram method (filter NOT applied):\t$tmp\" >> $fout.tmp";
		#print "Running command:\n$command\n\n";
		system($command);
	}
	system("cat $fout.tmp $fin > $fout");
	unlink "$fout.tmp";
}
else{
	if($filter_channel){
		$SNR_cutoff[2-$filter_channel] = 0;
	}
	my $SNRcutoffString = join(" ", @SNR_cutoff);
	$fout =~ s/\.bnx$//i;
	$command = "$RefAligner -i $fin -merge -bnx -minSNR $SNRcutoffString -o $fout -f\n";
	system($command);
}
exit 0;
######################################################
#                    Subroutines                     #
######################################################
sub Init{
	my $opt_string = 'hi:o:C:m:M:d:KL:P:';
	if(!getopts("$opt_string", \%opts)){
		print STDERR ("ERROR: Invalid parameter(s)! Try -h for more information.\n");
		Usage();
	}
	Usage() if $opts{h};
}

sub Check_inputs{
	my ($fin_filename, $fout_filename, $in_dir, $out_dir);
	if(!defined($opts{i})){
		print STDERR ("ERROR: Missing parameter (-i)! Try -h for more information.\n");
		Usage();
	}
	if( !($opts{i} =~ /\.bnx$/i) ){
		print STDERR ("ERROR: The input file must have a suffix of \"bnx\"!\n");
		exit (1);
	}
	else{
		if(! -f $opts{i}){
			print STDERR ("ERROR: The input file does not exist!\n");
			exit (1);
		}
		($fin_filename, $in_dir) = fileparse($opts{i});
		$in_dir = abs_path($in_dir);
		$fin = File::Spec->catfile($in_dir, $fin_filename);
		#print "IN: [$in_dir] [$fin_filename] [$fin]\n";
	}
	if(defined($opts{o})){
		($fout_filename, $out_dir) = fileparse($opts{o});
		$out_dir = abs_path($out_dir);
		$fout = File::Spec->catfile($out_dir, $fout_filename);
		#print "OUT: [$out_dir] [$fout_filename] [$fout]\n";
	}
	if(defined($opts{P})){
		#$fout_SNR = File::Spec->catfile($opts{P}, "SNR_histogram");
		$fout_SNR = $opts{P};
		$fout_SNR =~ s/\.pdf$//i;
	}
	if(defined($opts{C})){
		$filter_channel = $opts{C};
	}
	if(defined($opts{m})){
		$default_min_SNR = $opts{m};
	}
	if(defined($opts{M})){
		$default_max_SNR = $opts{M};
	}
	if(defined($opts{d})){
		$default_SNR = $opts{d};
	}
	if(defined($opts{L})){
		$min_num_of_labels = $opts{L};
	}
	#if(defined($opts{o})){
	#	$fout_summary = $fout;
	#	$fout_summary =~ s/\.bnx$/_summary.txt/i;
	#}
	#else{
		$OUT_summary = *STDOUT;
	#}

	my $file_is_good = `tail -20 $fin | grep '^QX' | wc -l`;
	chomp($file_is_good);
	if(! $file_is_good){
		#open($OUT_summary, ">".$fout_summary) || die ("ERROR: Can't open $fout_summary: $!\n");
		print $OUT_summary("ERROR! Unsupported file format!\n");
		#close($OUT_summary);
		exit 1;
	}
}

sub Usage{
	print << "EOF";

Usage: $^X $0 [Options] <Args>
Options:
  -h    : This help message
  -i    : Input BNX file (Required, can be 1 or 2 colors)
  -o    : Output BNX file (default: does not generate new BNX)
  -C    : Only filter the BNX on the channel specified (use with -o)
  -K    : Output BNX file BUT KEEP ALL labels (default: OFF, valid with -o)
  -m    : Minimum allowed dynamic SNR cutoff (default: 2)
  -M    : Maximum allowed dynamic SNR cutoff (default: 20)
  -d    : Default dynamic SNR cutoff (3.5) if outside the range of [m, M]
  -L    : Minimum number of labels for the script to estimate SNR
  -P    : Plot SNR histogram using the PDF filename specified (no default)
EOF
	exit (1);
}

sub GetSNR{
	my ($fin, $fout_SNR) = @_;
	
	my @hasQX;
	my $N_channels;
	my @SNR_array;
	my (@BIN, @BIN_less, @BIN_plus);
	my @step;
	my ($default_min_SNR_less, $default_max_SNR_plus) = ($default_min_SNR - 0.5, $default_max_SNR + 0.5);
	if($default_min_SNR_less <= 1){
		print STDERR ("ERROR: The input SNR cutoff is out of bound!\n");
		exit (1);
	}
	if($default_min_SNR >= $default_max_SNR){
		print STDERR ("ERROR: The Minimum SNR cutoff should be less than the Maximum SNR cutoff!\n");
		exit (1);
	}
	#my @min = ($default_max_SNR) x 2;
	my @max = ($default_min_SNR) x 2;
	my @total_labels = (0, 0);
	my (@SNR_count, @SNR_count_less, @SNR_count_plus);
	my $OUT_SNR;
	
	for(my $i = 0; $i < 2; $i++){
		my $QX_string = "QX". ($i+1) . "1";
		$hasQX[$i] = `tail -n 20 $fin | grep '^$QX_string' | head -n 1 | wc -l`;
		chomp($hasQX[$i]);
	}
	$N_channels = $hasQX[1] == 0 ? 1 : 2;
	for(my $i = 0; $i < $N_channels; $i++){
		$BIN[$i] = int( log($default_max_SNR / $default_min_SNR) * 25 );
		$step[$i] = log($default_max_SNR_plus / $default_min_SNR_less) / $BIN[$i];
	}
	
	open(my $IN, $fin) || die ("ERROR: Can't open $fin: $!\n");
	while(my $line = <$IN>){
		chomp $line;
		$line =~ s/\r//g;
		for(my $i = 0; $i < $N_channels; $i++){
			my $QX_string = "QX". ($i+1) . "1";
			if($line =~ /^$QX_string/){
				my @x = split("\t", $line);
				$total_labels[$i] += (@x-1);
				for(my $j = 1; $j < @x; $j++){
					#if($x[$j] < $min[$i]){
					#	$min[$i] = $x[$j];
					#}
					if($x[$j] > $max[$i]){
						$max[$i] = $x[$j];
					}
					# The middle range bins of the histogram
					if($x[$j] > $default_min_SNR_less && $x[$j] < $default_max_SNR_plus){
						#push(@{$SNR_array[$i]}, $x[$j]);
						my $idx = floor( log( $x[$j]/$default_min_SNR_less ) / $step[$i] );
						$SNR_count[$i]->[$idx]++;
					}
					# The left range bins of the histogram
					elsif($x[$j] <= $default_min_SNR_less && $x[$j] > 0){
						my $idx = ceil( log( $default_min_SNR_less/$x[$j] ) / $step[$i] );
						$SNR_count_less[$i]->[ceil( log($default_min_SNR_less) / $step[$i] )-$idx]++;
					}
					# The right range bins of the histogram
					elsif($x[$j] >= $default_max_SNR_plus){
						my $idx = floor( log( $x[$j]/$default_min_SNR_less ) / $step[$i] );
						$SNR_count_plus[$i]->[$idx-$BIN[$i]]++;
					}
				}
			}
		}
	}
	close($IN);
	
	if( defined($fout_SNR) ){
		open($OUT_SNR, ">$fout_SNR") || die ("ERROR: Can't open $fout_SNR: $!\n");
	}
	my @SNR_cutoff = ("NA", "NA");
	for(my $j = 0; $j < $N_channels; $j++){
		$BIN_less[$j] = ceil( log($default_min_SNR_less) / $step[$j] );
		$BIN_plus[$j] = ceil( log($max[$j]/$default_max_SNR_plus) / $step[$j] );
		if($BIN_less[$j] < 0){
			$BIN_less[$j] = 0;
		}
		if($BIN_plus[$j] < 0){
			$BIN_plus[$j] = 0;
		}
		
		if($hasQX[$j]){
			if($total_labels[$j] < $min_num_of_labels){
				$SNR_cutoff[$j] = -1;
			}
			else{
				#($SNR_cutoff[$j], $SNR_count[$j]) = getSNR_cutoff_from_SNRArray($SNR_array[$j], $min[$j], $max[$j], $BIN[$j]);
				for(my $i = 0; $i < $BIN[$j]; $i++){
					if( !defined($SNR_count[$j]->[$i]) ){
						$SNR_count[$j]->[$i] = 0;
					}
				}
				for(my $i = 0; $i < $BIN_less[$j]; $i++){
					if( !defined($SNR_count_less[$j]->[$i]) ){
						$SNR_count_less[$j]->[$i] = 0;
					}
				}
				for(my $i = 0; $i < $BIN_plus[$j]; $i++){
					if( !defined($SNR_count_plus[$j]->[$i]) ){
						$SNR_count_plus[$j]->[$i] = 0;
					}
				}
				
				my $local_min = 1e15;
				my $idx = 0;
				my @SNR_cutoff_tmp;
				my $end;
				for($end = $BIN[$j] - 1; $end >= 0; $end--){
					if( $SNR_count[$j]->[$end] == 0 ){
						next;
					}
					last;
				}
				my $head = 1;
				for(my $i = 0; $i <= $end; $i++){
					if( $head && $SNR_count[$j]->[$i] == 0){
						next;
					}
					else{
						$head = 0;
					}
					if($SNR_count[$j]->[$i] < $local_min){
						$local_min = $SNR_count[$j]->[$i];
						@SNR_cutoff_tmp = ( exp(log($default_min_SNR_less)+($i+0.5)*$step[$j]) );
						$idx = 1;
					}
					elsif($SNR_count[$j]->[$i] == $local_min){
						$SNR_cutoff_tmp[$idx++] = exp(log($default_min_SNR_less)+($i+0.5)*$step[$j]);
					}
				}
				$SNR_cutoff[$j] = percentile(\@SNR_cutoff_tmp, 75);
				
				# For plotting purpose
				if( defined($fout_SNR) ){
					print $OUT_SNR ("# Estimated SNR cutoff for channel ", $j+1, ":\t", $SNR_cutoff[$j], "\n");
					for(my $i = 0; $i < scalar(@{$SNR_count_less[$j]}); $i++){
						my $log_value = log($default_min_SNR_less) + ($i+0.5 - $BIN_less[$j]) * $step[$j];
						print $OUT_SNR ($j+1, "\t", log($default_min_SNR_less) + ($i - $BIN_less[$j]) * $step[$j], "\t", $log_value, "\t", exp($log_value), "\t", $SNR_count_less[$j]->[$i], "\n");
					}
					my $limit = scalar(@{$SNR_count[$j]});
					if(!defined($SNR_count_plus[$j])){
						$limit = $end + 1;
					}
					for(my $i = 0; $i < $limit; $i++){
						my $log_value = log($default_min_SNR_less) + ($i+0.5) * $step[$j];
						print $OUT_SNR ($j+1, "\t", log($default_min_SNR_less) + ($i) * $step[$j], "\t", $log_value, "\t", exp($log_value), "\t", $SNR_count[$j]->[$i], "\n");
					}
					if(defined($SNR_count_plus[$j])){
						for(my $i = 0; $i < scalar(@{$SNR_count_plus[$j]}); $i++){
							my $log_value = log($default_max_SNR_plus) + ($i+0.5) * $step[$j];
							print $OUT_SNR ($j+1, "\t", log($default_max_SNR_plus) + ($i) * $step[$j], "\t", $log_value, "\t", exp($log_value), "\t", $SNR_count_plus[$j]->[$i], "\n");
						}
					}
				}
			}
		}
	}
	if( defined($fout_SNR) ){
		close($OUT_SNR);
		
		my $script_dir = abs_path(dirname($0));
		my $command = "Rscript $script_dir/SNR_histogram_plot.R $fin $fout_SNR $fout_SNR.pdf " . join(" ", @SNR_cutoff);
		#print "$command\n";
		system($command);
	}
	return @SNR_cutoff;
}

sub getSNR_cutoff_from_SNRArray{
	my ($SNR_array_ref, $min, $max, $BIN) = @_;
	my @SNR_array = @$SNR_array_ref;
	my @SNR_count;
	my $step = (log($max) - log($min)) / $BIN;
	
	for(my $i = 0; $i < @SNR_array; $i++){
		my $idx = int( log($SNR_array[$i]/$min) / $step );
		$SNR_count[$idx]++;
	}
	splice @SNR_count, $BIN;  # We only use the first $BIN bins
	for(my $i = 0; $i < $BIN; $i++){
		if(! defined($SNR_count[$i]) ){
			$SNR_count[$i] = 0;
		}
	}
	
	my $local_min = 1e15;
	my $idx = 0;
	my @SNR_cutoff_tmp;
	for(my $i = 0; $i < $BIN; $i++){
		if($SNR_count[$i] < $local_min){
			$local_min = $SNR_count[$i];
			@SNR_cutoff_tmp = ( exp(log($min)+($i+0.5)*$step) );
			$idx = 1;
		}
		elsif($SNR_count[$i] == $local_min){
			$SNR_cutoff_tmp[$idx++] = exp(log($min)+($i+0.5)*$step);
		}
	}

	#for(my $i = 0; $i < @SNR_count; $i++){
	#	print exp(log($min)+($i+0.5)*$step), "\t", $SNR_count[$i], "\n";
	#}
	#my $result = percentile(\@SNR_cutoff_tmp, 75);
	#print "$result\n";
	return (percentile(\@SNR_cutoff_tmp, 75), \@SNR_count);
}

sub percentile{
	my ($array_ref, $percentile) = @_;
	my $count = @$array_ref;
	
	if(!$count){
		return undef;
	}
	if($percentile == 0){
		return $array_ref->[0];
	}
	if($percentile == 100){
		return $array_ref->[$count-1];
	}
	
	my $quantile_int = ($percentile/100) * ($count-1) + 1;  # 1-based
	my $quantile_float = $quantile_int - floor($quantile_int);
	$quantile_int = floor($quantile_int);
	
	# interpolation
	my $quantile_value = $array_ref->[$quantile_int-1];
	if($quantile_float == 0){
		return $quantile_value;
	}
	my $quantile_nextvalue = $array_ref->[$quantile_int];
	
	return ($quantile_value + $quantile_float * ($quantile_nextvalue - $quantile_value));
}

sub RunCommand{
	my ($command, $err_msg) = @_;
	$err_msg ||= "ERROR: Cannot excute the command:\n$command\n";
	if(system($command)){
		print STDERR ($err_msg);
	}
}

__END__



