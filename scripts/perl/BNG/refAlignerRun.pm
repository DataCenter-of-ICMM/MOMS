# $Id: refAlignerRun.pm 4044 2015-08-13 23:50:43Z apang $

package BNG::refAlignerRun;

use strict;
use Cwd qw(abs_path cwd);
use File::Basename;
use IO::Select;
use IPC::Open3;

our $VERSION = 1.00;
	
our @acceptedParams = ("binPath", "f", "if", "i", "ref", "o", "endoutlier", "outlier", "extend", "FN", "FP", "sf", "sd", "sr",
                       "res", "resSD", "mres", "A", "biaswt", "M", "Mfast", "maxmem", "maxthreads", "deltaX", "deltaY",
                       "xmapchim", "RepeatMask", "RepeatRec", "T", "BestRef", "nosplit", "XmapStatWrite", "merge", 
                       "stdout", "stderr", "first", "pairmerge", "preferNGS", "readparameters", "NoBpp", "minsites", "maxsites", "pairmergeRepeat");
our $switchParams = {
   	  		"f" 			=> 1,
   	  		"stdout" 		=> 1,
   	  		"stderr"		=> 1,
   	  		"merge" 		=> 1,
   	  		#"pairmerge"		=> 1, 
   	  		"preferNGS"		=> 1,
			"NoBpp"			=> 1,
			"pairmergeRepeat"	=> 1
   	    };
   	    
our $acceptedParams_href={};

## This adds "${CURRENT_SCRIPT_PATH}/perl5/" direcory at run time to the @INC array
## This script sould sit one level above the additional Perl modules directory "pe".
BEGIN {
  my $script_path = abs_path(dirname($0));
  my $module_path = abs_path($script_path . "/../../");
#  my $module_path_bno = abs_path(dirname($0) . "/../../pm");
#  # my $config_path = abs_path(dirname($0) . "/../hybridScaffold.pl");
  unshift @INC, $module_path;
}

{
    for (my $i=0; $i<=$#acceptedParams; $i++) {
   	  if (exists ($acceptedParams_href->{$acceptedParams[$i]})) { 
   	  	  print "\nWarning the refAligner option <$acceptedParams[$i]> is in the list already\n";
   	  } else {
   	  	 $acceptedParams_href -> {$acceptedParams[$i]} = 1;
   	  }
    }
}

sub new {
  my $class = shift;
  my $self = {};
  my $params = shift;
  $self->{binPath} = $params->{binPath};
  bless($self, $class);
  $self->setParams($params);
  return $self;
}

sub getBinPath { 
	my $self = shift;
	return $self->{binPath};
}

sub setBinPath {
	my $self = shift;
	$self->{binPath} = shift;
}

sub getParam { 
	my $self = shift;
	my $pp = shift;
   	if (exists ($self->{params}->{$pp}) ) { 
   	  	return (0, $self->{params}->{$pp});
   	} else { 
   	  	return (1, "Error: not defined parameter $pp");
   	}		
}

sub getParams { 
	my $self = shift;
    return $self->{params};		
}

sub rmParam { 
	my $self = shift;
	my $pp = shift;
   	if (exists ($self->{params}->{$pp}) ) { 
   	  	delete $self->{params}->{$pp};
   	} 	
}

sub setDefaultParams { 
	my $self = shift;
	my $dp = { 		
	};
	$self->setParams($dp);
}

sub setParams { 
   my $self = shift;
   my $userParams = shift;

   for my $p (keys %$userParams) {
   	  if (exists ($acceptedParams_href->{$p}) ) { 
   	  	if ($p ne 'binPath') {
   	  		$self->{params}->{$p} = $userParams->{$p};
   	    }
   	  } else { 
   	  	die "ERROR: not accepted refAligner parameter key <$p> with value as $userParams->{$p}\n";
   	  }
   }
}

sub runCMD {
	my $self = shift;
	my $outResults="";
	my $errResults="";
    my $cmd = $self->getCMD();
    # my $outResults = `$cmd`;
    # my $job_status=$? >> 8;
	my ($outResults, $errResults, $job_status) = $self->runCommand($self->getCMD_ARef());
	if ($job_status!=0) { 
	  $errResults = "ERROR: in calling $cmd with exit code $job_status and info: $outResults; " . $errResults;
	  warn "ERROR: $errResults\n";
	}	
	return ($outResults, $errResults, $job_status);
}

# get full command string: 
sub getCMD { 
	my $self = shift;
	my @pp = sort keys %{ $self->{params} };
	my $cmd = '"' . $self->getBinPath() . '"';
	for (my $i=0; $i <= $#pp; $i++) {
		my $pk = $pp[$i];
		if (exists $switchParams->{$pk}) {
			# print "$pk swith\n"; 
			$cmd = $cmd . " -" . $pk; 
		} else {
			my $val = $self->{params}->{$pk};
			if (ref $val) { 
				if ( ref $val eq "ARRAY" ) {
					# arrary value:
					my $nr = @$val;
					for (my $j=0; $j < $nr; $j++) { 
	            				$cmd = $cmd . " -" . $pk . " " . $val->[$j];
	            			}
				} else {
					my $ty = ref $val; 
					die "ERROR: Input parameter <$pk> has unacceptable data type $ty";
				}
			} else { 
              			# 'This is a scalar':
              			$cmd = $cmd . " -" . $pk . " " . $val;            	
            		} # if ref
		}
	}
	return $cmd;
}

# get full command string as an array reference: 
sub getCMD_ARef { 
	my $self = shift;
	my @pp = sort keys %{ $self->{params} };
	my $cmd = [$self->getBinPath()];
	for (my $i=0; $i <= $#pp; $i++) {
		my $pk = $pp[$i];
		if (exists $switchParams->{$pk}) {
			push(@$cmd, "-" . $pk); 
		} else {
			my $val = $self->{params}->{$pk};
			if (ref $val) { 
				if ( ref $val eq "ARRAY" ) {
					# arrary value:
					my $nr = @$val;
	            	for (my $j=0; $j < $nr; $j++) { 
	            		# push(@$cmd, " -" . $pk . " " . $val->[$j]);
	            		push(@$cmd, "-" . $pk);
	            		push(@$cmd, $val->[$j]);
	            	}
				} else {
					my $ty = ref $val; 
					die "ERROR: Input parameter <$pk> has unacceptable data type $ty";
				}
            } else { 
              # 'This is a scalar':
              #push(@$cmd, " -" . $pk . " " . $val); 
              push(@$cmd, "-". $pk);
              if ($pk eq "RepeatRec" || $pk eq "RepeatMask" || $pk eq "pairmerge") {
                 my (@tmp) = split(/\s+/, $val);
                 for (my $itmp = 0; $itmp<=$#tmp; $itmp++) {
                   push(@$cmd, $tmp[$itmp]);         
                 }  
	      		  } else { 
                 push(@$cmd, $val);
              }	
            }
		}
	}
	return $cmd;
}

sub runCommand  {
	    my $self = shift;
        my ($argsRef) = @_;
        # my $pwd = cwd();
        # my $tmp = join(", ", @$argsRef);
        # print "runCommand - arg: $tmp\n";
        my $pid = open3(my $CMD_IN, my $out, my $err, @$argsRef) || die "ERROR: Unable to execute ";
        # print "runCommand - pid=$pid\n";

        close($CMD_IN);

        my $outResults = "";
        my $errResults = "";
        my $sel = new IO::Select;
        $sel->add($out, $err);
        while(my @fhs = $sel->can_read) {
                foreach my $fh (@fhs) {
                        my $line = <$fh>;
                        unless(defined $line) {
                                $sel->remove($fh);
                                next;
                        } # unless line
                        if($fh == $out) {
                                $outResults .= "$line";
                                #print "$line";
                        }elsif($fh == $err) {
                                $errResults .= "$line";
                                #print "$line";
                        }else{
                                die "ERROR: This should never execute!";
                        } # if fh
                } # foreach fh
        } # while
        my $ret=waitpid ($pid, 0); # reap the exit code
        my $child_exit_status = $? >> 8;
        # print "runCommand - outResults=$outResults,\nerrResults=$errResults\n";
        return ($outResults, $errResults, $child_exit_status);
} # runCommand


1;

__END__
