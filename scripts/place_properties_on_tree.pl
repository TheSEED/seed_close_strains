########################################################################
# 
# This program takes as input a file describing the presence/absence of
# families in each of the genomes (the tips of the tree).  It produces
# the families.on.tree file, which extends the presence/absence data to
# nonleaf nodes.  The properties in that file are just
# 
# 	[Property,Node,Values]
# 
# where Values contains
# 
#      	     0   - does not include the family
# 	     1   - does include the family
# 	     0,1 - uncertain
#
use Carp;
use strict;
use Data::Dumper;
use Getopt::Long;
use SeedEnv;
use tree_utilities;
my $usage = "usage: place_properties_on_tree -t tree -p property-table -e extended-property-table\n";
my $treeF;
my $propF;
my $extended_propF;

my $rc  = GetOptions('t=s' => \$treeF,
		     'p=s' => \$propF,
		     'e=s' => \$extended_propF);

if ((! $rc) ||
    (! -s $treeF) || (! -s $propF) || (! open(EX,">$extended_propF")))
{ 
    print STDERR $usage; exit ;
}

my $tree;
($tree = &parse_newick_tree(join("",`cat $treeF`)))
    || die "could not parse the tree";
my $tips = &tips_of_tree($tree);
open(PROP,"<$propF") || die "could not open $propF";
my %prop_occ;
my %genomes;
while (defined($_ = <PROP>))
{
    chop;
    my($prop,$tip,$val) = split(/\t/,$_);
    $genomes{$tip} = 1;
    $prop_occ{$prop}->{$tip} = $val;
}
close(PROP);
my $tree = &subtree($tree,\%genomes);
my @props = keys(%prop_occ);
foreach my $prop (@props)
{
    my $occ_on_nodes = &place_prop_on_nodes($tree,$prop_occ{$prop});
    foreach my $node (keys(%$occ_on_nodes))
    {
	my $vals = $occ_on_nodes->{$node};
	my $set = join(",",@$vals);
	print EX join("\t",($prop,$node,$set)),"\n";
    }
}

close(EX);

sub place_prop_on_nodes {
    my($tree,$occ_on_tips) = @_;

    my $occH = {};
    foreach my $tip (keys(%$occ_on_tips))
    {
	my $val = defined($occ_on_tips->{$tip}) ? $occ_on_tips->{$tip} : 'unknown';	
	$occH->{$tip} = [$val];
    }
    &pass1($tree,$occH);
    &pass2($tree,$occH);
    return $occH;
}

sub pass1 {
    my($node,$occH) = @_;

    if (@{$node->[2]} > 1)  # if not child
    {
	my $i;
	my @current_values;
	for ($i=1; ($i < @{$node->[2]}); $i++)
	{
	    &pass1($node->[2]->[$i],$occH);
	    push(@current_values,$occH->{$node->[2]->[$i]->[0]});
	}
	$occH->{$node->[0]} = &iu(\@current_values);
    }
}

sub pass2 {
    my($node,$occH) = @_;

    if (@{$node->[2]} > 1)  # if not child
    {
	my $i;
	for ($i=1; ($i < @{$node->[2]}); $i++)
	{
	    my $overlap  = &i([$occH->{$node->[0]},$occH->{$node->[2]->[$i]->[0]}]);
	    if (@$overlap > 0)
	    {
		$occH->{$node->[2]->[$i]->[0]} = $overlap;
	    }
	    &pass2($node->[2]->[$i],$occH);
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
#    print &Dumper($values);
    foreach my $v (@$values)
    {
	if (! ref($v)) { confess('BAD') }
	foreach my $x (@$v)
	{
	    $counts{$x}++;
	}
    }
    my @tmp = sort grep { $counts{$_} == @$values } keys(%counts);
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
    return [sort keys(%counts)];
}

