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
use warnings;

# Show normal message
sub infoMsg
{
	my ($msg, $topic, $color, $strWidth) = @_;

	exit unless( defined($msg) );

	if( !defined $topic ){
		print "$msg"; # prefer no "\n"
		return;
	}

	if( $strWidth ){
		$topic = sprintf("%${strWidth}s", $topic);
	}

	if( $color ){
		print "\033[${color}m$topic: \033[0m$msg\n";
	}
	else{
		print "$topic: $msg\n";
	}
}

# Show error message
sub errMsg
{
	my ($msg, $topic) = @_;

	exit unless( defined($msg) );

	$topic = 'ERROR' unless(defined($topic));
	if( lc($topic) =~ /^warn/ ){
		print "\033[33mWarning: \033[0m$msg\n";
	}
	else{
		print "\033[31m$topic: \033[0m$msg\n";
	}
}

sub dieMsg
{
	my ($msg, $topic) = @_;

	$topic = 'ERROR' unless(defined($topic));
	die("\033[31m$topic: \033[0m$msg\n") if(defined($msg));
}

1;
