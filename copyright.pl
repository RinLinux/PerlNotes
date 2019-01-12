# 
# **/ 
use strict; 
use warnings; 

$^I = ".bak";
while (<>) {
	if (/\A# \//) {
		$_ .= "## Copyright (C) 2016 by Zhenglin\n";
	}
	print;
}


#perl copyright.pl email_regx.pl
