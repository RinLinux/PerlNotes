

use strict;
use warnings;

use File::Spec;
use Cwd;
use Data::Dumper;

my $curdir = getcwd();
foreach (sort glob(".* *")) {
	#print "    ", File::Spec->catfile($curdir,$_),"\n";
}

my $hash = {
	"a" => [1,2,3],
	"b" => (4,5,6),
	"c" => (7,8,9),
};

print $$hash{"a"}->[0],"\n";
print Dumper($hash);
