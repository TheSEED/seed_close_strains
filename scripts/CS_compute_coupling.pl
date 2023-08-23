use strict;
use Data::Dumper;
use Getopt::Long;
use SeedEnv;

my $usage = "usage: CS_compute_coupling -d DataDir\n";
my $dataD;
my $rc  = GetOptions('d=s' => \$dataD);

if ((! $rc) || (! -d $dataD))
{ 
    print STDERR $usage; exit ;
}

my @fam_func_peg = map { chomp; [split(/\t/,$_)] } `cut -f1,2,4 $dataD/families.all`;
my %fam_of_peg  = map { ($_->[2] => $_->[0]) } @fam_func_peg;
my %func_of_fam  = map { ($_->[0] => $_->[1]) } @fam_func_peg;
my @genomes = map { ($_ =~ /^(\S+)/) ? $1 : () } `cat $dataD/genome.names`;

my %close_fams;
foreach my $g (sort @genomes)
{
    print STDERR "processing $g\n";  ## handles only single reg
    my @locs = map { (($_ =~ /^(\S+)\t(\S+:)?(\S+)_(\d+)([+-])(\d+)$/) && $fam_of_peg{$1}) ? [$1,$3,$4,$5,$6] : () }
#               `echo $g | svr_all_features peg | svr_location_of`;
                `cat '$dataD/PegLocs/$g'`;
    @locs = sort { ($a->[1] cmp $b->[1]) or ($a->[2] <=> $b->[2]) } @locs;
    my $i;
    for ($i=0; ($i < $#locs); $i++)
    {
	my $j;
	for ($j = $i+1; ($j < @locs) && 
	                (($locs[$i]->[1] eq $locs[$j]->[1]) && 
                         (abs(($locs[$j]->[2] - $locs[$i]->[2])) < 6000)); $j++)
	{
	    my $f1 = $fam_of_peg{$locs[$i]->[0]};
	    my $f2 = $fam_of_peg{$locs[$j]->[0]};
	    if (defined($f1) && defined($f2))
	    {
		$close_fams{$f1}->{$f2}++;
		$close_fams{$f2}->{$f1}++;
	    }
	}
    }
}

open(COUPLED,">$dataD/coupled.families") || die "could not open $dataD/coupled.families";
foreach my $f1 (sort { $a <=> $b } keys(%close_fams))
{
    my $xH = $close_fams{$f1};
    foreach my $f2 (sort { $xH->{$b} <=> $xH->{$a} } keys(%$xH))
    {
	my $y = $xH->{$f2};
	if ($y > 1)
	{
	    print COUPLED join("\t",($f1,$f2,$y)),"\n";
	}
    }
}
close(COUPLED);
