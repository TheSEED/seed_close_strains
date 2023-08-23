use Carp;
use strict;
use Data::Dumper;
use Getopt::Long;
use SeedEnv;
use tree_utilities;
my $usage = "usage: where_shifts_occurred -t tree -e extended-property-table\n";
my $treeF;
my $extended_propF;

my $rc  = GetOptions('t=s' => \$treeF,
		     'e=s' => \$extended_propF);

if ((! $rc) ||
    (! -s $treeF) || (! -s $extended_propF))
{ 
    print STDERR $usage; exit ;
}

my $tree;
($tree = &parse_newick_tree(join("",`cat $treeF`)))
    || die "could not parse the tree";

open(PROP,"<$extended_propF") || die "could not open $extended_propF";
my %properties;
while (defined($_ = <PROP>))
{
    chop;
    my($prop,$node,$val) = split(/\t/,$_);
    $properties{$node}->{$prop} = $val;
}
close(PROP);
&show_shifts($tree,\%properties);

sub show_shifts {
    my($node,$props) = @_;

    if (@{$node->[2]} > 1)  # if not child
    {
	my $props_for_node = $properties{$node->[0]};
	my $i;
	for ($i=1; ($i < @{$node->[2]}); $i++)
	{
	    my $props_for_desc = $properties{$node->[2]->[$i]->[0]};
	    foreach my $p (keys(%$props_for_node))
	    {
		my $parent_prop = $props_for_node->{$p};
		my $desc_prop   = $props_for_desc->{$p};
		if ($parent_prop ne $desc_prop)
		{
		    print join("\t",($node->[0],$node->[2]->[$i]->[0],$p,$parent_prop,$desc_prop)),"\n";
		}
	    }
	    &show_shifts($node->[2]->[$i],$props);
	}
    }
}

	    
