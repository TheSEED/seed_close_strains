########################################################################
use strict;
use Data::Dumper;
use Getopt::Long;
use SeedEnv;

#
# This program produces a file of adjacency "properties".
# 
# 338	+	2000	downstream	1140.7	fig|1140.7.peg.193	fig|1140.7.peg.192
# 
# is a typical output line.  It says that family 338 has downstream adjacency
# to plus (same) strand instance of family 2000 in genome 1140.7 due to
# pegs fig|1140.7.peg.193 and fig|1140.7.peg.192.  For each family that occurs just
# once in the genome, there will be two such "properties", assuming the adjacent family
# occurs just once and the midpoints of the two genes differ by less than 5000 bp.
#
# NOTE: This construction of adjacency just ignores genes in sets in which a single
#       genome has multiple members.  This leads to "skipped genes" that look
#       perfectly conserved, while the difference is within 5kb but a ways off.
#        
#       I have not yet formed an opinion of whether or not this is the
#       behavior we want.
#        
#        
my $usage = "usage: cs_adjacency_data -d DataDir\n";
my $dataD;
my $rc  = GetOptions('d=s' => \$dataD);
if ((! $rc) || (! -d $dataD))
{ 
    print STDERR $usage; exit ;
}

my %peg_to_fam;
my %peg_to_func;
open(FAM,"<$dataD/families.all") || die "could not open $dataD/families.all";
my $last = <FAM>;
while ($last && ($last =~ /^(\S+)/))
{
    my $fam = $1;
    my @set;
    my %genomes;
    my $mult = 0;
    while ($last && ($last =~ /^(\S+)\t([^\t]*)\t[^\t]*\t(\S+)/) && ($1 eq $fam))
    {
	my $func = $2;
	my $peg = $3;
	push(@set,[$peg,$func]);
	my $g = &SeedUtils::genome_of($peg);
	if ($genomes{$g})
	{
	    $mult=1;
	}
	else
	{
	    $genomes{$g} = 1;
	}
	$last = <FAM>;
    }
    if (! $mult)
    {
	foreach my $tuple (@set)
	{
	    my($peg,$func) = @$tuple;
	    $peg_to_fam{$peg} = $fam;
	    $peg_to_func{$peg} = $func;
	}
    }

}


opendir(LOCS,"$dataD/PegLocs") || die "could not open $dataD/PegLocs";
my @genomes = grep { $_ !~ /^\./ } readdir(LOCS);
closedir(LOCS);

open(OUT,"| sort -T . -k 1n -k 4 > $dataD/adjacency.of.unique") || die "could not open $dataD/adjacency.of.unique";
foreach my $g (@genomes)
{
    my @tuples;
    foreach $_ (`cat '$dataD/PegLocs/$g'`)
    {
	if (($_ =~ /^(\S+)\t([^:]*:)?(\S+)_(\d+)([-+])(\d+)$/) && $peg_to_fam{$1})
	{
	    my($peg,$contig,$b,$strand,$ln) = ($1,$3,$4,$5,$6);
	    my $e = ($strand eq '+') ? ($b + ($ln-1)) : ($b - ($ln-1));
	    push(@tuples,[$peg,[$contig,$b,$e],$peg_to_fam{$peg},$peg_to_func{$peg}]);
	}
    }
    @tuples = sort { ($a->[1]->[0] cmp $b->[1]->[0]) or
		     (($a->[1]->[1] + $a->[1]->[2]) <=> ($b->[1]->[1] + $b->[1]->[2])) } @tuples;
    my $i;
    for ($i=0; ($i < @tuples); $i++)
    {
	my($c1,$b1,$e1) = @{$tuples[$i]->[1]};
	if ($i > 0)
	{
	    my($c2,$b2,$e2) = @{$tuples[$i-1]->[1]};
	    if (((($b1+$e1) - ($b2+$e2))/2) <= 5000)
	    {
		if ($b1 < $e1)
		{
		    my $dir2 = ($b2 < $e2) ? '+' : '-';
		    print OUT join("\t",($tuples[$i]->[2],
					 $dir2,
					 $tuples[$i-1]->[2],
					 'upstream',
					 $g,
					 $tuples[$i]->[0],
					 $tuples[$i-1]->[0])),"\n";
		}
		else
		{
		    my $dir2 = ($b2 < $e2) ? '-' : '+';
		    print OUT join("\t",($tuples[$i]->[2],
					 $dir2,
					 $tuples[$i-1]->[2],
					 'downstream',
					 $g,
					 $tuples[$i]->[0],
					 $tuples[$i-1]->[0])),"\n";
		}
	    }
	}
	if ($i < $#tuples)
	{
	    my($c2,$b2,$e2) = @{$tuples[$i+1]->[1]};
	    if (((($b2+$e2) - ($b1+$e1))/2) <= 5000)
	    {
		if ($b1 < $e1)
		{
		    my $dir2 = ($b2 < $e2) ? '+' : '-';
		    print OUT join("\t",($tuples[$i]->[2],
					 $dir2,
					 $tuples[$i+1]->[2],
					 'downstream',
					 $g,
					 $tuples[$i]->[0],
					 $tuples[$i+1]->[0])),"\n";
		}
		else
		{
		    my $dir2 = ($b2 < $e2) ? '-' : '+';
		    print OUT join("\t",($tuples[$i]->[2],
					 $dir2,
					 $tuples[$i+1]->[2],
					 'upstream',
					 $g,
					 $tuples[$i]->[0],
					 $tuples[$i+1]->[0])),"\n";
		}
	    }
	}
    }
}
close(OUT);

open(IN,"<","$dataD/adjacency.of.unique") || die "could not open $dataD/adjacency.of.unique";
open(OUT,">","$dataD/events") || die "could not open events";

my $last = <IN>;
while ($last && ($last =~ /^(\d+)\t[+-]\t\d+\t(\S+)/))
{
    my $fam = $1;
    my $dir1 = $2;
    my @set;
    while ($last && ($last =~ /^(\d+)\t[+-]\t\d+\t(\S+)/) && ($1 == $fam) && ($2 eq $dir1))
    {
	push(@set,$last);
	$last = <IN>;
    }
    &process(\@set,\*OUT,\%peg_to_func);
}
close(OUT);
close(IN);

use tree_utilities;
(my $rooted = &parse_newick_tree(join("",`cat $dataD/labeled.tree`)))
    || die "could not parse the tree";
my ($traits,$node_to_traits) = &place_events_on_leaves($dataD);
foreach my $trait (@$traits)
{
    &place_trait_on_tree($trait,$node_to_traits,$rooted);
}

open(OUT,">$dataD/traits.on.nodes") || die "could not open $dataD/traits.on.nodes";
foreach my $node (sort keys(%$node_to_traits))
{
    my $traitH = $node_to_traits->{$node};
    foreach my $trait (sort keys(%$traitH))
    {
	my $v = $node_to_traits->{$node}->{$trait};
	if (! $v) { die "HERE" }
	print OUT join("\t",($node,$trait,join(",",sort(@$v)))),"\n";
    }
}
close(OUT);

open(OUT,"| sort -T . > $dataD/placed.events") || die "could not open $dataD/placed.events";
&placed_events($rooted,$node_to_traits,\*OUT);
close(OUT);

sub placed_events {
    my($tree,$node_to_traits,$fh) = @_;

    my $this_node = $tree->[0];
    my $this_node_traits = $node_to_traits->{$this_node};
    my @traits_for_this_node = grep { @{$this_node_traits->{$_}} == 1 } keys(%$this_node_traits);
    my $i;
    for ($i=1; ($i < @{$tree->[2]}); $i++)
    {
	&placed_events($tree->[2]->[$i],$node_to_traits,$fh);
	my $child = $tree->[2]->[$i]->[0];
	my $child_traits = $node_to_traits->{$child};
	foreach my $trait (@traits_for_this_node)
	{
	    my $child_val;
	    if (($child_val = $child_traits->{$trait}) && (@$child_val == 1))
	    {
		my $this_node_val = $this_node_traits->{$trait}->[0];
		if ($child_val->[0] ne $this_node_val)
		{
		    print $fh join("\t",($this_node,$child,$trait,$this_node_val,$child_val->[0])),"\n";
		}
	    }
	}
    }
}

sub place_trait_on_tree {
    my($trait,$node_to_trait,$tree) = @_;
    my $values = {};
    &pass1($trait,$tree,$node_to_trait,$values);
    &pass2($trait,$tree,$node_to_trait,$values);
    foreach my $node (keys(%$values))
    {
	if (@{$values->{$node}} > 0) 
	{
	    $node_to_trait->{$node}->{$trait} = $values->{$node};
	}
    }
}

sub pass1 {
    my($trait,$tree,$node_to_trait,$values) = @_;

    if (@{$tree->[2]} == 1)
    {
	my $v = $node_to_trait->{$tree->[0]}->{$trait};
	if (! $v) { $v = [] }
	$values->{$tree->[0]} = $v;
    }
    else
    {
	my $i;
	my @current_values;
	for ($i=1; ($i < @{$tree->[2]}); $i++)
	{
	    &pass1($trait,$tree->[2]->[$i],$node_to_trait,$values);
	    push(@current_values,$values->{$tree->[2]->[$i]->[0]});
	}
	$values->{$tree->[0]} = &iu(\@current_values);
    }
}

sub pass2 {
    my($trait,$tree,$node_to_trait,$values) = @_;

    if (@{$tree->[2]} == 1)
    {
	my $v = $node_to_trait->{$tree->[0]}->{$trait};
	if (! $v) { $v = [] }
	$values->{$tree->[0]} = $v;
    }
    else
    {
	my $i;
	for ($i=1; ($i < @{$tree->[2]}); $i++)
	{
	    my $overlap  = &i([$values->{$tree->[0]},$values->{$tree->[2]->[$i]->[0]}]);
	    if (@$overlap > 0)
	    {
		$values->{$tree->[2]->[$i]->[0]} = $overlap;
	    }
	    &pass2($trait,$tree->[2]->[$i],$node_to_trait,$values);
	}
    }
}

sub iu {
    my($values) = @_;

    my $x = &i($values);
    if (@$x > 0)
    {
	return $x;
    }
    return &u($values);
}

sub i {
    my($values) = @_;
    my %counts;
    foreach my $v (@$values)
    {
	if (! ref($v)) { confess('BAD') }
	foreach my $x (@$v)
	{
	    $counts{$x}++;
	}
    }
    my @tmp = grep { $counts{$_} == @$values } keys(%counts);
    return \@tmp;
}
	    
sub u {
    my($values) = @_;

    my %counts;
    foreach my $v (@$values)
    {
	foreach my $x (@$v)
	{
	    $counts{$x}++;
	}
    }
    return [keys(%counts)];
}

sub place_events_on_leaves {
    my($dataD) = @_;

    my $node_to_traits = {};
    my %traits;
    foreach my $tuple ( map { my($f1,$o2,$f2,$o1,$g) = split(/\t/,$_); [$f1,$o1,$f2,$o2,$g] } `cat $dataD/events`)
    {
	my($f1,$o1,$f2,$o2,$g) = @$tuple;
	my $trait = join(":",($f1,$o1));
	my $val   = join(":",($f2,$o2));
	$node_to_traits->{$g}->{$trait} = [$val];
	$traits{$trait} = 1;
    }
    my $traits = [sort keys(%traits)];
    return($traits,$node_to_traits);
}

sub process {
    my($set,$fh,$to_func) = @_;

    $set->[0] =~ /^\d+\t([+-]\t\d+\t\S+)/;
    my $key = $1;
    my $i;
    for ($i=1; ($i < @$set) && &same($key,$set->[$i]); $i++) {}
    if ($i < @$set)
    {
	foreach $_ (@$set)
	{
	    chop;
	    my @flds = split(/\t/,$_);
	    push(@flds,$to_func->{$flds[5]},$to_func->{$flds[6]});
	    print $fh join("\t",@flds),"\n";
	}
    }
}

sub same {
    my($key,$line) = @_;
    $line =~ /^\d+\t([+-]\t\d+\t\S+)/;
    return ($1 eq $key);
}

