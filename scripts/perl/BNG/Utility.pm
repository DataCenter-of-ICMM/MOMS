package BNG::Utility;

use strict;
use Cwd qw(abs_path);
use POSIX qw[strftime ceil];
use List::Util qw[min max];

BEGIN{
	use Exporter ();
	use vars qw ($AUTHOR $VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
	$AUTHOR	= 'Institute of Chinese Materia Medica, China Academy of Chinese Medical Sciences';
	@EXPORT	= qw(
			readCMap
			writeCMapFile
			getContigStat
			getSubsetCMap
			shiftCMapIds
			countOverhangLabels
			readXMap
			writeXMapFile
			parsingFastMrgStdOut
			writeAllFastMrgPairs
			parseConfig
			makeParams
			alignDP
			alignLCS
			alignGappedDP
			compareScore
			locate
			calculateLinearParams
			getPerfectLinearParams
				);
	@EXPORT_OK	 = qw();
	%EXPORT_TAGS = ();
	@ISA 		 = qw(Exporter);
	$VERSION	 = 0.01;
}

##
# read a CMAP file
##
sub readCMap
{
	my ($cmap_file) = @_;
	# read cmap file (first time to see if there is the keyword "ChimQuality")
	open(IN, "$cmap_file") or die ("ERROR: Unable to read in file $cmap_file: $!\n");
	my $foundChimQuality = 0;
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/)	{
			# header line
			if ($line =~ /#h\s+/){
				$foundChimQuality = 1 if ($line =~ /ChimQuality/);
				last;	# done
			}
		} # if line
	} # while line
	close IN;
	# then read it again, will read in according whether ChimQuality is there or not
	
	my $cmap={};
	my $c_cmapId = 0;
	my $numcontig = 0;
	my @contigLength = ();
	$cmap->{"FileName"} = abs_path($cmap_file);
	my @header = ();
	my @data_name = ();
	my @data_type;
	my $cmap_version = "";
	my $nchannels;
	my %channels = ();
	open(IN, "$cmap_file") or die ("ERROR: Unable to read in file $cmap_file: $!\n");
	my $NUM_C_LIMITED = ($foundChimQuality == 1) ? (10) : (9);	# assumption: column 10 has ChimQuality, if it is present
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		if ($line =~ /^#/)	{
			# header line
			if ($line =~ /^#h\s+/)	{
				# data name:
				my $headTag = "";	# the #h
				($headTag, @data_name) = split(/\s+/, $line);
				splice @data_name, $NUM_C_LIMITED;
				push(@header, join("\t", "$headTag $data_name[0]", @data_name[1..$#data_name]));
				next;
			}
			if ($line =~ /^#f\s+/)	{
				# data type:
				my $headTag = "";       # the #f
				($headTag, @data_type) = split(/\s+/, $line);
				splice @data_type, $NUM_C_LIMITED;
				push(@header, join("\t", "$headTag $data_type[0]", @data_type[1..$#data_type]));
				next;
			}
			if ($line =~ /^#\s+CMAP File Version:\s+(\S+)/)	{
				$cmap_version = $1;
			} elsif($line =~ /^#\s+Label Channels:\s+(\S+)/){
				$nchannels = $1;
			} elsif($line =~ /^#\s+Nickase Recognition Site (\d+):\s+(\S+)/){
				$channels{$1} = $2;
			} 
			push(@header, $line);
			next;
		} # if header line
		if($nchannels != scalar(keys %channels)){
			die("ERROR: inconsistent channel number found in the header of $cmap_file\n");
		}
		
		# data lines
		if (! exists $cmap->{"headers"})	{
			# first data line, populate the relevant hash values
			$cmap->{"headers"} = \@header;
			$cmap->{"dataName"} = \@data_name;
			$cmap->{"dataType"} = \@data_type;
			$cmap->{"channels"} = \%channels;
			$cmap->{"version"} = $cmap_version;
		}
		
		my $numc = $NUM_C_LIMITED;
		
		# now store information of that data line
		$line =~ s/^\s*//;	$line =~ s/\s+/\t/g;
		my @d_items = split(/\t/, $line);
		if ($d_items[0] != $c_cmapId)	{
			# a new contig:
			$numcontig += 1;
			$c_cmapId=$d_items[0];
			# ContigLength:
			$cmap->{contigs}->{$c_cmapId}->{$data_name[1]} = $d_items[1];
			push(@contigLength, $d_items[1]);
			# NumSites:
			$cmap->{contigs}->{$c_cmapId}->{$data_name[2]} = $d_items[2];
			
			# the rest of the columns
			for (my $i = 3; $i < $numc; $i += 1)	{
				$cmap->{contigs}->{$c_cmapId}->{$data_name[$i]} = [$d_items[$i]];
			} # for i
		} else	{
			# another site of the same contig:
			for (my $i = 3; $i < $numc; $i += 1)	{
				my $ad = $cmap->{contigs}->{$c_cmapId}->{$data_name[$i]};
				push(@$ad, $d_items[$i]);
			} # for i
		} # if same cmap id
	} # while line
	close IN;
	$cmap->{nContigs} = $numcontig;
	return ($cmap, $numcontig, \@contigLength);
}

sub getLargestCMapId
{
	my ($hash) = @_;
	my $big = 0;
	foreach(keys %$hash){
		if($_ > $big){
			$big = $_;
		}
	}
	return $big;
}

sub writeCMapFile 
{
	my ($cmap, $cmap_file, $want_full_header, $min_sites, $max_sites) = @_;

	$want_full_header = 1 unless(defined $want_full_header);

	#print "   write a cmap file $cmap_file ...\n";
	
	open CMF, ">$cmap_file" or die ("ERROR: Unable to write to file " + $cmap_file + " - " + $!);	
	
	my @data_name = @{ $cmap->{dataName} };
	my @cmapIds = sort {$a <=> $b} keys %{ $cmap->{contigs} };
	if(defined $min_sites or defined $max_sites){
		my $k = 0;
		my $no = scalar(@cmapIds);
		for (my $i=0; $i < $no; $i++) {
			my $c = $cmapIds[$i];
			my $numsites = $cmap->{contigs}->{$c}->{$data_name[2]};
			next if(defined $min_sites && $numsites < $min_sites);
			next if(defined $max_sites && $numsites > $max_sites);
			$cmapIds[$k++] = $c;
		}
		while($k<$no){
			pop(@cmapIds);
			$k++;
		}
	}
	## headers:
	my @headers = @{ $cmap->{headers} };
	if ($want_full_header) {
		for (my $i=0; $i<=$#headers; $i++) {
			if($headers[$i] =~ /^# (Number of Consensus .*[Mm]aps)/){
				printf CMF "# $1:\t%d\n", scalar(@cmapIds);
			} else{
				print CMF $headers[$i] . "\n";
			}
		}
	} else { 
		# only header and format:
		print CMF "$headers[-2]\n";
		print CMF "$headers[-1]\n";
	}
	for (my $i=0; $i < scalar(@cmapIds); $i++) {
		my $c = $cmapIds[$i];
		my $numsites = $cmap->{contigs}->{$c}->{$data_name[2]};
		for (my $j=0; $j <= $numsites; $j++) { 
			print CMF $c . "\t" . $cmap->{contigs}->{$c}->{$data_name[1]} . "\t" . $numsites;
			for (my $m=3; $m <= $#data_name; $m++) {
				print CMF  "\t" . $cmap->{contigs}->{$c}->{$data_name[$m]}->[$j];
			}
			print CMF "\n";
		}
	}
	
	close CMF;
}

##
# create a new cmap with some cmap ids are excluded:
##
sub getSubsetCMap
{
	my ($cmap, $givenCMapID_href, $flag) = @_;
	$flag = "exclude" unless(defined $flag);
	# also deal with list as an array ref
	if (ref $givenCMapID_href eq "ARRAY") { 
		my %givenCMapID = map { $_ => 1 } @$givenCMapID_href;
		$givenCMapID_href = \%givenCMapID;
	}
    my $allEID = join(",", keys %{$givenCMapID_href});	
	my $sCMap={};
	# - no header - $sCMap->{headers}
	print "   getting subset of cmap data from $cmap->{FileName}...\n";
    $sCMap->{dataName} = $cmap->{dataName};
    $sCMap->{dataType} = $cmap->{dataType};
    $sCMap->{version} = $cmap->{version};

    my @comment=("# This CMAP file was subset of $cmap->{FileName} with CMapID " . $flag .": " . $allEID);
    push(@comment, @{$cmap->{headers}});
    $sCMap->{headers} = \@comment;
#contigs:
    my $numc = 0;
    for my $contig ( sort {$a <=> $b} keys %{ $cmap->{contigs} })   { 
    	#print "contig - $contig ";
    	if ($flag =~ /excl/) { 
    		# exclude given ids:
	    	if ($givenCMapID_href->{$contig}) { 
	    		# print "$contig excluded *******\n";
	    	} else { 
	    		#print "included\n";
	    		$sCMap->{contigs}->{$contig} = $cmap->{contigs}->{$contig};
	    		$numc++;
	    	}
    	} else { 
    		# include given ids only:
	    	if ($givenCMapID_href->{$contig}) {
	    	  # print "included\n";
	    		$sCMap->{contigs}->{$contig} = $cmap->{contigs}->{$contig};
	    		$numc++; 
	    	} else { 
	    		# print "$contig excluded *******\n";
	    	}    		
    	}
    }
    $sCMap->{nContigs} = $numc;  
    return $sCMap;
}

sub getContigStat
{
 # give a contig length array, caculate its statistical data:
    my $len = @_;
    if ($len == 0) {
    	return ("NA", "NA", "NA", "NA", "NA", "NA");
    }
    my @vals = sort {$b <=> $a} @_;
    my $median;
    my $r = $len%2;
    my $m = ($len-$r)/2;
    if($r != 0) #odd?
    {
        $median = $vals[$m];
        #print "\n$r, $m\n";
    }
    else #even
    {
        $median = ($vals[$m-1] + $vals[$m])/2;
        #print "\n$r, $m\n";
    }
    my $min = $vals[-1];
    my $max = $vals[0];
    my $total = 0.0;
    for (my $i=0; $i<$len; $i++) { 
    	$total = $total + $vals[$i];
    }
    my $mean = $total/($len);

#	N50 length is $n50: 
#   N50 value is $n50value: 
#   L50 is $L50:
    my $n50=0; 
    my $L50=0;
    my $n50value=0.0;
    if ($len==1) {
    	$n50value = $max;
    	$n50=$max;
    	$L50=1;
    } else {
	    for (my $i=0; $i<$len; $i++) {
	     $n50=$n50+$vals[$i];
	     $L50++;
	     if($n50 == $total/2.0){
	     	$n50value=$vals[$i];
	     	last; 
	     }
	     if($n50 > $total/2.0){
	     	$n50value=($vals[$i-1]+$vals[$i])/2.0;
	        last; 
	     }
	    }
    }
   	# print "N50 length is $n50 and N50 value is: $n50value and L50 is $L50\n"; 
    return ($min, $max, $mean, $median, $n50value, $total);
}

sub shiftCMapIds
{
	my ($cmap, $nshift) = @_; 

	my @cmap_ids = sort {$b <=> $a} keys %{$cmap->{contigs}};
	for (my $i = 0; $i <= $#cmap_ids; $i++) { 
		my $contig_k = $cmap->{contigs}->{$cmap_ids[$i]};
		#delete the hash element:
		delete $cmap->{contigs}->{$cmap_ids[$i]};
		$cmap->{contigs}->{ $cmap_ids[$i] + $nshift } = $contig_k;
	}
}
##
# read xmap file:
##
sub readXMap
{
	my ($xmap_file) = @_;
	my $in;
	if($xmap_file eq "-"){
		$in = \*STDIN;
	}
	else{
		open XMF, "<$xmap_file" or die ("ERROR: Unable to read in file " + $xmap_file + " - " + $!);
		$in = \*XMF;
	}
	my $xmap={};
	$xmap->{"FileName"} = abs_path($xmap_file);
	my @header=();
	my @data_name;
	my @data_type;
	my $xmap_version;
	my $nchannels;
	#processing header/comment lines:
	my $line=<$in>;
	while ($line =~ /^#/) {
		chomp($line);
		# bug fix in case sometime col-name is RefcontigID in some earlier verisons (should be RefContigID):
		if ($line =~ /^#h\s+/) {
			if ($line =~ /RefcontigID/) { 
				$line =~ s/RefcontigID/RefContigID/;
			}
		}
		push (@header, $line);
		if ($line =~ /^#h\s+/) { 
			# data name:
			(undef, @data_name) = split(/\s+/, $line);
		}
		if ($line =~ /^#f\s+/) { 
			# data type:
			(undef, @data_type) = split(/\s+/, $line);
			last;
		}
		if ($line =~ /^#\s+XMAP File Version:/) {
			my (@tmp) = split(/\s+/, $line);
			$xmap_version = $tmp[-1];
		}
		elsif($line =~ /^#\s+Label Channels:\s+(\S+)/){
			$nchannels = $1;
		}
		$line=<$in>;		
	}
	$xmap->{"headers"} = \@header;
	$xmap->{"dataName"} = \@data_name;
	$xmap->{"dataType"} = \@data_type;
	$xmap->{"version"} = $xmap_version;
	$xmap->{"nchannels"} = $nchannels;

	#processing data:
	#XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum	QryLen RefLen	LabelChannel	Alignment
	my $numc = scalar(@data_name);
	for (my $i=0; $i<$numc; $i++) { 
		$xmap->{hits}->{$data_name[$i]} = [];
	}
	$line=<$in>;
	my $c_xmapId=0;
    my $numm = 0;
	while (defined $line) {
		if($line =~ /^#/){
			$line = <$in>;
			next;
		}
		# print $line;
		chomp($line);
		my (@d_items) = split(/\s+/, $line);
		$numm++;
        for (my $i=0; $i<$numc; $i++) {
        	my $a = $xmap->{hits}->{$data_name[$i]};
        	push(@$a, $d_items[$i]);
        	$xmap->{hits}->{$data_name[$i]} = $a;
        }
		$line=<$in>;
	}
	$xmap->{totalHits} = $numm;
	close $in;

	return ($xmap);	
}

##
# write a xmap data to a file
##
sub writeXMapFile
{
	my ($xmap, $xmap_file, $wantHeader, $rCmapFileName, $qCmapFileName) = @_;

	$wantHeader = 1 unless(defined $wantHeader);

	my $out;
	if(defined $xmap_file){
		open XMF, ">$xmap_file" or die ("ERROR: Unable to write in file " + $xmap_file + " - " + $!);
		$out = \*XMF;
	}
	else{
		$out = \*STDOUT;
	}

	if ($wantHeader == 1)	{
		# print the header lines
		my $tmp = $xmap->{headers};
		my $numLines = @$tmp;
		for (my $i = 0; $i < $numLines; $i++)	{
			my $line = $xmap->{headers}->[$i];
			if ($line =~ /^#\s+Reference/ && defined $rCmapFileName)	{
				print $out "# Reference Maps From:\t$rCmapFileName\n";
			} elsif ($line =~ /^#\s+Query/ && defined $qCmapFileName)	{
				print $out "# Query Maps From:\t$qCmapFileName\n";
			} else	{
				print $out "$line\n";
			} # if line
		} # for i
	} # if wantHeader
	
   	my @data_name = @{ $xmap->{dataName} };
   	my $numc = @data_name;
	for (my $i=0; $i < $xmap->{totalHits}; $i++) { 
		for (my $j=0; $j < $numc; $j++) { 
	  		print $out $xmap->{hits}->{$data_name[$j]}->[$i];
	  		if ($j < $numc - 1) { 
	  			print $out "\t";
	  		} else { 
	  			print $out "\n";
	  		}
	  	}
	}
	close $out;
}

##
# count number of labels extending outside a region:
##
sub countOverhangLabels
{
	my ($cmap, $cmapId, $leftPosition, $rightPosition) = @_;

	my ($numL, $numR) = (0, 0);
	my $Position = $cmap->{contigs}->{$cmapId}->{Position};
	my $totalNS = $cmap->{contigs}->{$cmapId}->{NumSites};
	for (my $i=0; $i<$totalNS; $i++) { 
		if ( $Position->[$i] < $leftPosition ) { 
			$numL++;
		};
		if ( $Position->[$i] > $rightPosition ) { 
			$numR++;
		};
	}
	return ($numL, $numR);
}

sub parsingFastMrgStdOut
{
	my ($stdout_file) = @_;
	open (IN, $stdout_file) or die ("ERROR: Unable to read in file " + $stdout_file + " - " + $!);
	my $merge_pairs = {};
	my @contigId1 = ();
	my @contigId2 = ();
	my @resultsContigId = ();
	my @mergeRounds = ();
	my $numMrgP = 0;
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		my ($theStage, $firstId, $firstSize, $secondId, $secondSize) = (-1, -1, -1, -1, -1);
		my ($hybridId, $hybridSize) = (-1, -1);
		if ($line =~ /^(\d+).+Creat.* merged map \(id=(\d+),len=(\d+\.\d*)\) based on alignment between .*\(id=(\d+),len=(\d+\.\d*)\).*\(id=(\d+),len=(\d+\.\d*)\)/)	{
			# parsing out these lines Created merged map (id=785,len=200.000) based on alignment between Y mapid=784(id=785,len=200.0) and X mapid=42(id=43,len=100.000):mapid=784,id=785,len=200.0,sites=20 (Just keeping Ymap=0xffffffff)
			# or
			# Created merged map (id=43, len=250.000) based on alignment between Y mapid=784(id=785,len=200.0) and X mapid=42(id=43,len=100.000):mapid=1001,id=43,len=250.0,sites=25
			($theStage, $hybridId, $hybridSize, $firstId, $firstSize, $secondId, $secondSize) = ($1, $2, $3, $4, $5, $6, $7);

			push(@mergeRounds, $theStage);
			push(@contigId1, $firstId);
			push(@contigId2, $secondId);
			push(@resultsContigId, $hybridId);
			$numMrgP += 1;

		} elsif ($line =~ /^(\d+).+Creat.* merged map based on alignment between .*\(id=(\d+),len=(\d+\.\d*)\).*\(id=(\d+),len=(\d+\.\d*)\).+/)	{
			# parsing out from these lines (\d+):Created merged map based on alignment between .*\(id=(\d+),len=.*\(id=(\d+),len=.*\:id=(\d+),len=.*
			($theStage, $firstId, $firstSize, $secondId, $secondSize) = ($1, $2, $3, $4, $5);
			if ($line =~ /\(Just keeping Xmap.+\)$/ || $line =~ /\(Just keeping Ymap.+\)$/)       {
				# one of the maps is completely encompassed by another map, then the larger map's id is the final hybrid id
				$hybridId = ($firstSize < $secondSize) ? ($secondId) : ($firstId);
			} else	{
				# hybrid Id is the smaller of the two ids
				$hybridId = ($firstId < $secondId) ? ($firstId) : ($secondId);
			} # if line has that string
			push(@mergeRounds, $theStage);
			push(@contigId1, $firstId);
			push(@contigId2, $secondId);
			push(@resultsContigId, $hybridId);
			$numMrgP += 1;

		} else	{
			# no op
		}
	}
	close IN;
	$merge_pairs->{ContigID1} = \@contigId1;
	$merge_pairs->{ContigID2} = \@contigId2;
	$merge_pairs->{ResultContigID} = \@resultsContigId;
	$merge_pairs->{numMrgPairs} = $numMrgP;
	$merge_pairs->{mergeRounds} = \@mergeRounds;

	return ($numMrgP, $merge_pairs);
}

sub writeAllFastMrgPairs
{
	my ($fn, $mrg_pairs, $numMrgP) = @_;

	my @mrg_rounds_ids = ("A".."Z", "AA".."AZ", "BA".."BZ", "CA".."CZ", "DA".."DZ");
	my $mrgRoundsIdsRef = \@mrg_rounds_ids;
	open(OUT, ">$fn") or die("ERROR: Unable to write to file " + $fn + " - " + $!);
	print OUT '"RowID" "ContigID1" "ContigID2" "ResultContigID" "MrgRound"'."\n";
	for (my $i = 0; $i < $numMrgP; $i++)	{
		my $letterMrgRound = ($mrg_pairs->{mergeRounds}[$i] >= scalar(@$mrgRoundsIdsRef)) ? ($mrg_pairs->{mergeRounds}[$i]) : ($mrgRoundsIdsRef->[$mrg_pairs->{mergeRounds}[$i] - 1]);	
		print OUT '"'.($i + 1).'" "'.$mrg_pairs->{ContigID1}[$i].'" "'.
			$mrg_pairs->{ContigID2}[$i].'" "'.
			$mrg_pairs->{ResultContigID}[$i].'" "Mrg'.
			$letterMrgRound."\"\n";	
	}
	close OUT;
}

sub parseConfig
{
	my ($configRef, $stage) = @_;
	my %paraDict = ();
	my $cnt = 0;
	my @configs = ();
# only support the 'include' of level one
	my $config = $configRef->{$stage};
	foreach my $element (keys %{$config}){
		if($element =~ /include/){
			my $eletype = ref($config->{$element});
			my $attributes;
			my @Attributes = ();
			if($eletype eq "HASH"){
				push(@Attributes, $config->{$element});
				$attributes = \@Attributes;
			}
			elsif($eletype eq "ARRAY"){
				$attributes = \@Attributes;
			}
			for(my $i=0; $i<scalar(@{$attributes}); $i++){
				foreach my $val (keys %{$attributes->[$i]}){
					if($val =~ /val(\d+)/){
						push(@configs, $configRef->{$attributes->[$i]{$val}});
					}
				}
			}
		}
	}
	push(@configs, $configRef->{$stage});
	foreach $config (@configs){
		foreach my $element (keys %{$config}){
			my $eletype = ref($config->{$element});
			if($element =~ /flag/){
				my $attributes;
				my @Attributes = ();
				if($eletype eq "HASH"){
					push(@Attributes, $config->{$element});
					$attributes = \@Attributes;
				}
				elsif($eletype eq "ARRAY"){
					$attributes = \@{$config->{$element}}
				}
				for(my $i=0; $i<scalar(@{$attributes}); $i++){
					my $tmpAttr;
					my %tmpVals = ();
					foreach my $flagkey (keys %{$$attributes[$i]}){
						if($flagkey =~ /attr/){
							$tmpAttr = $$attributes[$i]{$flagkey};
							$tmpAttr =~ s/^-//; 
						}
						elsif($flagkey =~ /val(\d+)/){
							$tmpVals{$1} = $$attributes[$i]{$flagkey};
						}
					}
					if(defined $tmpAttr){
						my @tmpArr = ();
						foreach my $idx (sort {$a <=> $b} keys %tmpVals){
							push(@tmpArr, $tmpVals{$idx});
						}
						$paraDict{$tmpAttr}{val} = join(" ", @tmpArr);
						$paraDict{$tmpAttr}{idx} = $cnt++;
					}
				}
				next;
			}
		}
	}

	return \%paraDict;
}

sub makeParams
{
	my ($map) = @_;
	my @array = ();
	foreach my $attr ( sort { $map->{$a}{idx} <=> $map->{$b}{idx} } keys %{$map} ){
		push(@array, "-$attr");
		push(@array, "$map->{$attr}{val}");
	}
	return join(" ", @array);
}

sub alignDP
{
	my ($positions1, $positions2, $bidirectional) = @_;

	my $delta = 3; # number of successive label mismatches allowed
	my $alpha = 0.05; # expected variation ratio of fragment lengths
# initialization
	my @scores = (); # three dimentional score including number of matched fragments, matching scores, and the number of indel
	my ($i, $j, $k, $l);
	my ($m, $n) = ($#{$positions1}, $#{$positions2});
	my $maxIndel = max(ceil(min($m, $n) * $alpha), (defined $bidirectional ? 3 : 2));
	$scores[0][0] = [0, 0, 0];
	if(defined $bidirectional){
		$j = 1;
		my $iRowVal = ($bidirectional & 0x1) ? 0 : 9**9;
		# first row
		while($j <= $n - $m){
			$scores[0][$j++] = [0, 0, 0];
		}
		while($j <= $n){
			$scores[0][$j++] = [0, $iRowVal, 0];
		}
		my $iColVal= ($bidirectional & 0x2) ? 0 : 9**9;
		# first column
		$i = 1;
		while($i <= $m - $n){
			$scores[$i++][0] = [0, 0, 0];
		}
		while($i <= $m){
			$scores[$i++][0] = [0, $iColVal, 0];
		}
	}
	else{
		# first row
		$j = 1;
		while($j <= $n - $m){
			 $scores[0][$j++] = [0, 0, 0];
		}
		while($j <= $n){
			$scores[0][$j++] = [0, 9**9, 0];
		}
		# first column
		for($i=1; $i<=$m; $i++){
			$scores[$i][0] = [0, 0, 0];
		}
	}
	my @backPointer = ();
	my ($q, $r);
	my ($dev, $nl, $nid, $sc);
	my $score;
# fill the score matrix
	for($i=1; $i<=$m; $i++){
		for($j=1; $j<=$n; $j++){
			$scores[$i][$j] = [0, 9**9, 0];
			# for each cell, search for its (delta+1)^2 ajacent cells for an optimal path
			for($k=max($i-$delta-1, 0); $k<=$i-1; $k++){
				$q = $positions1->[$i] - $positions1->[$k];
				for($l=max($j-$delta-1, 0); $l<=$j-1; $l++){
					$score = $scores[$k][$l];
					$nid = ($i - $k - 1) + ($j - $l - 1);
					next if($score->[2] + $nid > $maxIndel);
					$r = $positions2->[$j] - $positions2->[$l];
					$dev = abs(($q - $r) / ($q + $r));
					$nl = (($dev <= $alpha) ? ($score->[0] + 1) : $score->[0]);
					$sc = $score->[1] + ($dev / $alpha) ** 2 + $delta * $nid;
					if(compareScore([$nl, $sc], $scores[$i][$j]) > 0){
						$scores[$i][$j] = [$nl, $sc, $score->[2] + $nid];
						$backPointer[$i][$j] = [$k, $l];
					}
				}
			}
		}
	}
	my $maxScore = $scores[$m][$n];
	my $pointer = [$m, $n];
	for($j=$n-1; $j>=0; $j--){
		if(compareScore($maxScore, $scores[$m][$j]) < 0){
			$maxScore = $scores[$m][$j];
			$pointer = [$m, $j];
		}
	}
	if( ($m > $n) or ($bidirectional & 0x01) ){
		for($i=$m-1; $i>=0; $i--){
			if(compareScore($maxScore, $scores[$i][$n]) < 0){
				$maxScore = $scores[$i][$n];
				$pointer = [$i, $n];
			}
		}
	}
	my @alignment = ();
	my $score = 0;
	if($maxScore->[0] > 0){
		my $prevPointer = $backPointer[$pointer->[0]][$pointer->[1]];
		my @subAlignments = ();
		my $status = 0;
		while(defined $prevPointer){
			if($scores[$pointer->[0]][$pointer->[1]]->[0] > $scores[$prevPointer->[0]][$prevPointer->[1]]->[0]){
				if($status == 0){
					unshift(@subAlignments, [$pointer]);
				}
				unshift(@{$subAlignments[0]}, $prevPointer);
				$status++;
			}
			else{
				$status = 0;
			}
			$pointer = $prevPointer;
			$prevPointer = $backPointer[$pointer->[0]][$pointer->[1]];
		}
# find the optimal sub-alignment first
		my @sorted = sort{ scalar(@{$subAlignments[$b]}) <=> scalar(@{$subAlignments[$a]}) || $scores[$subAlignments[$a]->[-1]->[0]][$subAlignments[$a]->[-1]->[1]]->[1] - $scores[$subAlignments[$a]->[0]->[0]][$subAlignments[$a]->[0]->[1]]->[1] <=> $scores[$subAlignments[$b]->[-1]->[0]][$subAlignments[$b]->[-1]->[1]]->[1] - $scores[$subAlignments[$b]->[0]->[0]][$subAlignments[$b]->[0]->[1]]->[1] } 0..$#subAlignments;
		@alignment = @{$subAlignments[$sorted[0]]};

# then extend the optimal sub-alignment
		my $k = $sorted[0] - 1;
		while($k >= 0){
			my $subAlign = $subAlignments[$k];
			my $safeGuardPointer = $subAlign->[-1];
			$pointer = $alignment[0];
			last if($scores[$pointer->[0]][$pointer->[1]]->[1] - $scores[$safeGuardPointer->[0]][$safeGuardPointer->[1]]->[1] > 16);
			$prevPointer = $backPointer[$pointer->[0]][$pointer->[1]];
			while($prevPointer->[0] != $safeGuardPointer->[0] or $prevPointer->[1] != $safeGuardPointer->[1]){
				unshift(@alignment, $prevPointer);
				$prevPointer = $backPointer[$prevPointer->[0]][$prevPointer->[1]];
			}
			unshift(@alignment, @{$subAlign});
			$k--;
		}
		$k = $sorted[0] + 1;
		while($k <= $#subAlignments){
			my @subAlign = @{$subAlignments[$k]};
			my $safeGuardPointer = $alignment[-1];
			$pointer = $subAlign[0];
			last if($scores[$pointer->[0]][$pointer->[1]]->[1] - $scores[$safeGuardPointer->[0]][$safeGuardPointer->[1]]->[1] > 16);
			$prevPointer = $backPointer[$pointer->[0]][$pointer->[1]];
			while($prevPointer->[0] != $safeGuardPointer->[0] or $prevPointer->[1] != $safeGuardPointer->[1]){
				unshift(@subAlign, $prevPointer);
				$prevPointer = $backPointer[$prevPointer->[0]][$prevPointer->[1]];
			}
			push(@alignment, @subAlign);
			$k++;
		}
		my $missed;
		if(defined $bidirectional){
			if(scalar(@alignment) < 5){ # to reduce false positive alignments
				return ([], 0);
			}
			if($bidirectional == 1){
				$missed = ($alignment[-1]->[1] == $n) ? $alignment[0]->[0] : ($alignment[0]->[0] + ($m - $alignment[-1]->[0]));
			}
			elsif($bidirectional == 2){
				$missed = ($alignment[0]->[1] == 0) ? ($m - $alignment[-1]->[0]) : ($alignment[0]->[0] + ($m - $alignment[-1]->[0]));
			}
			else{
				$missed = ($alignment[0]->[1] == 0) ? ($m - $alignment[-1]->[0]) :
						  ($alignment[-1]->[1] == $n) ? $alignment[0]->[0] : ($alignment[0]->[0] + ($m - $alignment[-1]->[0]));
			}
			if( $missed > scalar(@alignment) * 0.333 ){
				return ([], 0);
			}
		}
		else{
			$missed = ($m - $alignment[-1]->[0]) + $alignment[0]->[1];
			if( $missed >= (scalar(@alignment) * 2) ){
				return ([], 0);
			}
		}
		$score =  $scores[$alignment[-1]->[0]][$alignment[-1]->[1]]->[1] - $scores[$alignment[0]->[0]][$alignment[0]->[1]]->[1];
	}
	return (\@alignment, $score);
}

sub alignLCS
{
	my ($positions1, $positions2) = @_;

	my $delta = 3; # number of successive label mismatches allowed
	my $alpha = 0.05; # maximum allowed deviation dominated by total length for a pair of fragments to be regarded as matching
# initialization
	my @scores = (); # two dimentional score including number of matched labels and matching scores
	my ($i, $j, $k, $l);
	my ($m, $n) = ($#{$positions1}, $#{$positions2});
	for($j=0; $j<=$n; $j++){ # first row
		$scores[0][$j] = [0, 0];
	}
	for($i=1; $i<=$m; $i++){ # first column
		$scores[$i][0] = [0, 0];
	}
	my @backPointer = ();
	my ($q, $r);
	my ($dev, $nl, $sc);
	my $score;
# fill the score matrix
	for($i=1; $i<=$m; $i++){
		for($j=1; $j<=$n; $j++){
			$scores[$i][$j] = [0, 9**9];
			# for each cell, search for its (delta+1)^2 ajacent cells for an optimal path
			for($k=max($i-$delta-1, 0); $k<=$i-1; $k++){
				$q = $positions1->[$i] - $positions1->[$k];
				for($l=max($j-$delta-1, 0); $l<=$j-1; $l++){
					$r = $positions2->[$j] - $positions2->[$l];
					$dev = abs(($q - $r) / ($q + $r));
					$score = $scores[$k][$l];
					$nl = (($dev <= $alpha) ? ($score->[0] + 1) : $score->[0]);
					$sc = $score->[1] + ($dev / $alpha) ** 2 + $delta * (($i - $k - 1) + ($j - $l - 1));
					if($sc < $scores[$i][$j]->[1]){
						$scores[$i][$j] = [$nl, $sc];
						$backPointer[$i][$j] = [$k, $l];
					}
				}
			}
		}
	}
	my @alignment = ();
	my $k = 1;
	my $maxScore = $scores[$m][$k];
	for($j=2; $j<=$n; $j++){
		if(compareScore($maxScore, $scores[$m][$j]) < 0){
			$k = $j;
			$maxScore = $scores[$m][$k];
		}
	}
	my $pointer = [$m, $k];
	my $prevPointer = $backPointer[$pointer->[0]][$pointer->[1]];
	my $status = 0;
	while(defined $prevPointer){
		if($scores[$pointer->[0]][$pointer->[1]]->[0] > $scores[$prevPointer->[0]][$prevPointer->[1]]->[0]){
			unshift(@alignment, $pointer) if($status == 0);
			unshift(@alignment, $prevPointer);
			$status++;
		}
		else{
			$status = 0;
		}
		$pointer = $prevPointer;
		$prevPointer = $backPointer[$pointer->[0]][$pointer->[1]];
	}
	return \@alignment;
}

sub alignGappedDP
{
	my ($positions1, $positions2, $mi) = @_;

	my $delta = 3; # number of successive label mismatches allowed
	my $alpha = 0.05; # expected variation ratio of fragment lengths
	my $theta = 0.25; # expected variation ratio of gap size
# initialization
	my @scores = (); # three dimentional score including number of matched fragments, matching scores, and the number of indel
	my ($i, $j, $k, $l);
	my ($m, $n1, $n) = ($#{$positions1}, $mi-1, $#{$positions2});
	my $n2 = $n - $mi;
	my $maxIndel = ceil(min($m, $n1) * $alpha); # POSITIVELY related to SENSITIVIEY
	$scores[0][0] = [0, 0, 0];
	# first column
	$i = 1;
	while($i <= $m - $n1){
		$scores[$i++][0] = [0, 0, 0];
	}
	while($i <= $m){
		$scores[$i++][0] = [0, 9**9, 0];
	}
	my @rowIndex = (0);
	my @rowScore = (0);
	for($i=1; $i<=$m; $i++){
		my $fragSize  = $positions1->[$i] - $positions1->[$i-1];
		$rowIndex[$i] = $rowIndex[$i-1] + countScore($fragSize);
		$rowScore[$i] = $rowScore[$i-1] + sizeScore($fragSize);
	}
	my @colIndex = (0);
	my @colScore = (0);
	for($j=1; $j<=$n; $j++){
		my $fragSize  = $positions2->[$j] - $positions2->[$j-1];
		$colIndex[$j] = $colIndex[$j-1] + countScore($fragSize);
		$colScore[$j] = $colScore[$j-1] + sizeScore($fragSize);
	}
	my @backPointer = ();
	my ($q, $r);
	my ($dev, $nl, $nid, $sc, $rnid, $cnid);
	my $score;
	my $pointer;
	my ($minI, $minJ, $default) = (1, 0, 0);
# fill the score matrix
	for($j=1; $j<=$n; $j++){
		if($j == $mi){
			$r = $positions2->[$mi] - $positions2->[$mi-1];
			$minI = locate($positions1, $positions1->[0] + $r * (1 - $theta));
			$minI = 1 if($minI < 1);
			$scores[$minI-1][$j] = [0, 9**9, 0];
			for($i=$minI; $i<=$m; $i++){
				$scores[$i][$mi] = [0, 9**9, 0];
				my ($pos1, $pos2) = ($positions1->[$i] - $r * (1 + $theta), $positions1->[$i] - $r * (1 - $theta));
				my $iStart = locate($positions1, $pos1);
				my $iEnd = locate($positions1, $pos2);
				for($k=$iStart; $k<=$iEnd; $k++){
					$score = $scores[$k][$n1];
					$q = $positions1->[$i] - $positions1->[$k];
					$dev = abs(($q - $r) / ($q + $r));
					$nl = $score->[0];
					$sc = $score->[1] + ($dev / $alpha) ** 2;
					if(specialCompareScore([$nl, $sc], $scores[$i][$mi]) > 0){
						$scores[$i][$mi] = [$nl, $sc, 0];
						$backPointer[$i][$mi] = [$k, $n1];
					}
				}
			}
			$default = 9**9;
			$minJ = $mi;
			$maxIndel = ceil(min($m, $n2) * $alpha); # POSITIVELY related to SENSITIVITY
			next;
		}
		else{
			$scores[$minI-1][$j] = [0, $default, 0];
		}
		for($i=$minI; $i<=$m; $i++){
			$scores[$i][$j] = [0, 9**9, 0];
			# for each cell, search for its (delta+1)^2 ajacent cells for an optimal path
			for($k=max($i-$delta-1, $minI-1); $k<=$i-1; $k++){
				$q = $positions1->[$i] - $positions1->[$k];
				for($l=max($j-$delta-1, $minJ); $l<=$j-1; $l++){
					$score = $scores[$k][$l];
					$rnid = (($rowIndex[$i] > $rowIndex[$k]) ? (($rowIndex[$i] - $rowIndex[$k] - 1) * max($rowScore[$i] - $rowScore[$k] - 1, 0)): 0) ;
					$cnid = (($colIndex[$j] > $colIndex[$l]) ? (($colIndex[$j] - $colIndex[$l] - 1) * max($colScore[$j] - $colScore[$l] - 1, 0)): 0);
					$nid = $rnid + $cnid;
					next if($score->[2] + $nid > $maxIndel);
					$r = $positions2->[$j] - $positions2->[$l];
					$dev = abs(($q - $r) / ($q + $r));
					$nl = (($dev <= $alpha) ? ($score->[0] + 1) : $score->[0]);
					$sc = $score->[1] + ($dev / $alpha) ** 2 + $delta * $nid;
					if(compareScore([$nl, $sc], $scores[$i][$j]) > 0){
						$scores[$i][$j] = [$nl, $sc, $score->[2] + $nid];
						$backPointer[$i][$j] = [$k, $l];
					}
				}
			}
		}
	}
	my $maxScore = $scores[$m][$n];
	$pointer = [$m, $n];
	for($j=$n-1; $j>=$mi; $j--){
		if(compareScore($maxScore, $scores[$m][$j]) < 0){
			$maxScore = $scores[$m][$j];
			$pointer = [$m, $j];
		}
	}
	if($m > $n2){
		$minI = $n2 if($n2 > $minI);
		for($i=$m-1; $i>=$minI; $i--){
			if(compareScore($maxScore, $scores[$i][$n]) < 0){
				$maxScore = $scores[$i][$n];
				$pointer = [$i, $n];
			}
		}
	}
	my @alignment = ();
	my $score = [0, 9**9];
	if($maxScore->[0] > 0){
		my $lastPointer = $pointer;
		my $prevPointer = $backPointer[$pointer->[0]][$pointer->[1]];
		my @subAlignments = ();
		my $status = 0;
		while(defined $prevPointer){
			if( $scores[$pointer->[0]][$pointer->[1]]->[0] > $scores[$prevPointer->[0]][$prevPointer->[1]]->[0] or
				$pointer->[1] == $mi and ($scores[$pointer->[0]][$pointer->[1]]->[0] == 0 or $pointer->[0] == $m) ){
				if($status == 0){
					unshift(@subAlignments, [$pointer]);
				}
				unshift(@{$subAlignments[0]}, $prevPointer);
				$status++;
			}
			else{
				$status = 0;
			}
			$pointer = $prevPointer;
			$prevPointer = $backPointer[$pointer->[0]][$pointer->[1]];
		}
		my $firstPointer = $pointer;
# find the optimal sub-alignment first
		my @sorted = sort{ scalar(@{$subAlignments[$b]}) <=> scalar(@{$subAlignments[$a]}) || $scores[$subAlignments[$a]->[-1]->[0]][$subAlignments[$a]->[-1]->[1]]->[1] - $scores[$subAlignments[$a]->[0]->[0]][$subAlignments[$a]->[0]->[1]]->[1] <=> $scores[$subAlignments[$b]->[-1]->[0]][$subAlignments[$b]->[-1]->[1]]->[1] - $scores[$subAlignments[$b]->[0]->[0]][$subAlignments[$b]->[0]->[1]]->[1] } 0..$#subAlignments;
		@alignment = @{$subAlignments[$sorted[0]]};
		my ($hits1, $hits2);
		$hits1 = getHits1(\@alignment, $mi);
		$hits2 = getHits2(\@alignment, $mi);
		my %used = ($sorted[0] => 1);
		my $k;
		my $bSingleSide = 0;
		if($hits1 >= $hits2){
			if($hits2 == 0){
				$k = $sorted[0] + 1;
				while($k <= $#subAlignments){
					my $subAlign = $subAlignments[$k];
					$hits2 = getHits2($subAlign, $mi);
					if($hits2 > 0){
						if( ($k+1 <= $#subAlignments) and getHits2($subAlignments[$k+1], $mi) > $hits2 ){
							$k++;
							$subAlign = $subAlignments[$k];
						}
						push(@alignment, @{$subAlign});
						$used{$k} = 1;
						last;
					}
					$k++;
				}
				if($k > $#subAlignments){
					$bSingleSide = 1;
				}
			}
		}
		else{ # $hits2 > $hits1
			if($hits1 == 0){
				$k = $sorted[0] - 1;
				while($k >= 0){
					my $subAlign = $subAlignments[$k];
					$hits1 = getHits1($subAlign, $mi);
					if($hits1 > 0){
						if( ($k-1 >= 0) and getHits1($subAlignments[$k-1], $mi) > $hits1 ){
							$k--;
							$subAlign = $subAlignments[$k];
						}
						unshift(@alignment, @{$subAlign});
						$used{$k} = 1;
						last;
					}
					$k--;
				}
				if($k < 0){
					$bSingleSide = 1;
				}
			}
		}
		if($bSingleSide){
			return ([], [0, 9**9]); # obsolete it at present
		}

# then extend the optimal sub-alignment
		my @kk = (sort{$a <=> $b} keys %used);
		$k = $kk[0] - 1;
		while($k >= 0){
			my $subAlign = $subAlignments[$k];
			my $safeGuardPointer = $subAlign->[-1];
			$pointer = $alignment[0];
			last if($scores[$pointer->[0]][$pointer->[1]]->[1] - $scores[$safeGuardPointer->[0]][$safeGuardPointer->[1]]->[1] > 16);
			$prevPointer = $backPointer[$pointer->[0]][$pointer->[1]];
			while(defined $prevPointer and ($prevPointer->[0] != $safeGuardPointer->[0] or $prevPointer->[1] != $safeGuardPointer->[1])){
				unshift(@alignment, $prevPointer);
				$prevPointer = $backPointer[$prevPointer->[0]][$prevPointer->[1]];
			}
			unshift(@alignment, @{$subAlign});
			$k--;
		}
		$k = $kk[-1] + 1;
		while($k <= $#subAlignments){
			my @subAlign = @{$subAlignments[$k]};
			my $safeGuardPointer = $alignment[-1];
			$pointer = $subAlign[0];
			last if($scores[$pointer->[0]][$pointer->[1]]->[1] - $scores[$safeGuardPointer->[0]][$safeGuardPointer->[1]]->[1] > 16);
			$prevPointer = $backPointer[$pointer->[0]][$pointer->[1]];
			while(defined $prevPointer and ($prevPointer->[0] != $safeGuardPointer->[0] or $prevPointer->[1] != $safeGuardPointer->[1])){
				unshift(@subAlign, $prevPointer);
				$prevPointer = $backPointer[$prevPointer->[0]][$prevPointer->[1]];
			}
			push(@alignment, @subAlign);
			$k++;
		}
		if( $maxScore->[0] < max(0.7 * ($lastPointer->[1] - $firstPointer->[1]), 5) ){ # to reduce false positive hits
			return ([], [0, 9**9]);
		}
		my $final_score = $maxScore->[1]; # $scores[$alignment[-1]->[0]][$alignment[-1]->[1]]->[1] - $scores[$alignment[0]->[0]][$alignment[0]->[1]]->[1];
		if( $final_score > getScoreThreshold(scalar(@alignment)) ){ # to reduce false positive hits
			return ([], [0, 9**9]);
		}
		$score = [scalar(@alignment), $final_score];
	}
	return (\@alignment, $score);
}

sub compareScore
{
	my ($score1, $score2) = @_;
	if($score1->[0] < $score2->[0]){
		return -1;
	}
	if($score1->[0] > $score2->[0]){
		return 1;
	}
	if($score1->[1] > $score2->[1]){
		return -1;
	}
	if($score1->[1] < $score2->[1]){
		return 1;
	}
	return 0;
}

sub specialCompareScore
{
	my ($score1, $score2) = @_;
	my $diff = abs($score1->[0] - $score2->[0]);
	if($diff > 3){
		return ($score1->[0] < $score2->[0]) ? -1 : 1;
	}
	if($score1->[1] > $score2->[1] + 0.5 + $diff * 16){
		return -1;
	}
	if($score1->[1] < $score2->[1] - 0.5 + $diff * 16){
		return 1;
	}
	return compareScore($score1, $score2);
}

sub getHits1
{
	my ($align, $mi) = @_;

	return ($align->[0]->[1] < $mi) ? ( ($align->[-1]->[1] < $mi) ? scalar(@{$align}) : $mi - $align->[0]->[1] ) : 0;
}

sub getHits2
{
	my ($align, $mi) = @_;

	return ($align->[-1]->[1] > $mi) ? ( ($align->[0]->[1] > $mi) ? scalar(@{$align}) : $align->[-1]->[1] - $mi ) : 0;
}

sub sizeScore
{
	my ($fragSize) = @_;
	if($fragSize <= 600){
		return 0;
	}
	if($fragSize >= 2200){
		return 1;
	}
	return sprintf("%.4f", $fragSize * 0.000625 - 0.5);
}

sub getScoreThreshold
{
	my ($nAligned) = @_;
	if($nAligned <= 6){
		return 3 * $nAligned;
	}
	return sprintf("%.2f", ($nAligned * 1.05 - 3.3) * $nAligned);
}

sub countScore
{
	my ($fragSize) = @_;
	return ($fragSize > 600) ? 1 : 0;
}

sub locate
{
	my ($positions, $pos) = @_;
	my ($left, $right) = (0, $#{$positions});
	my ($middle, $val);
	while($left <= $right){
		$middle = int(($left + $right) / 2);
		$val = $positions->[$middle];
		if($pos < $val - 0.1){
			$right = $middle - 1;
			next;
		}
		if($pos > $val + 0.1){
			$left = $middle + 1;
			next;
		}
		return $middle;
	}
	return ($left == $middle) ? $middle : $left;
}

sub calculateLinearParams
{
	my ($xstart, $xend, $ystart, $yend) = @_;

	my $slope = ($ystart - $yend) / ($xstart - $xend);
	my $intercept = ($xstart * $yend - $xend * $ystart) / ($xstart - $xend);

	return ($slope, $intercept);
}

sub getPerfectLinearParams
{
	my ($xstart, $xend, $ystart, $yend) = @_;

	my $slope = (($xstart - $xend) * ($ystart - $yend) >= 0) ? 1 : -1;
	my $intercept = (($ystart + $yend) - $slope * ($xstart + $xend)) / 2;

	return ($slope, $intercept);
}

1;
