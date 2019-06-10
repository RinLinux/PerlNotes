

use Data::Dumper;

open my $coconet,"coconet.dat" or die $!;
my %totalbytes;
while (<$coconet>) {
	my ($source, $destination, $bytes) = split;
	$totalbytes{$source}{$destination} += $bytes;
}

# print Dumper(%totalbytes);
for my $source (sort keys %totalbytes){
	for my $des (sort keys %{$totalbytes{$source}}){
		print "$source => $des: \t$totalbytes{$source}{$des}\n";
	}
	print "\n";
}