use strict;
use Data::Dumper;

my $last = <STDIN>;
while ($last && ($last =~ /^PROTEIN-ID\s+(\S+)/))
{
    my $peg = $1;
    my @kmers;
    $last = <STDIN>;
    while ($last && ($last !~ /^PROTEIN-ID/))
    {
	if ($last =~ /^HIT\s+\d+\s+(\S+)/)
	{
	    push(@kmers,$1);
	}
	$last = <STDIN>;
    }
    print "$peg\t",join(",",sort { $a <=> $b } @kmers),"\n";
}
