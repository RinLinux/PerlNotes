


use File::Find;

sub get_dir {
	print "$File::Find::name was found \n";
}

find(\&get_dir,'/');

my $totall_size;

find(sub{ $totall_size += -s if -f },'.');
# print $totall_size;