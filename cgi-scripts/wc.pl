# -*- perl -*-
#
# Copyright (c) 2003-2006 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
# 
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License. 
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#

eval {
    require FIG_Config;
};

use URI::Escape;  # uri_escape
use gjoseqlib;
use HTML;
use strict;
use CGI;
our $cgi = new CGI;
use SeedEnv;
use tree_utilities;
use CloseStrains;

if (0) {
    print $cgi->header;
    my @params = $cgi->param;
    print "<pre>\n";
    foreach $_ (@params) {
		print "$_\t:",join(",",$cgi->param($_)),":\n";
    }
    print "\n",join('&',map { "$_=" . $cgi->param($_) } @params),"\n";
    exit;
}

our $html     = [];

our $node     = $cgi->param('node');
our $family   = $cgi->param('family');
our $genome   = $cgi->param('genome');
our $ali      = $cgi->param('alignment');
our $tree     = $cgi->param('tree');
our $request  = $cgi->param('request');
our $keywords = $cgi->param('keywords');
our $function = $cgi->param('function');
our $reaction = $cgi->param('reaction');
our $dataD    = $cgi->param('dataD');
our $base     = $cgi->param('base');
#
# If we are operating in RAST, the csD sits within the
# RAST job directory and the dataD setting is the
# RAST job number.
#

our $csD;
if ($base)
{
    $csD = $base;
}
else
{
    $csD      = "/homes/overbeek/Ross/MakeCS.Kbase/Data/CS";
}
if ($FIG_Config::rast_jobs && -d $FIG_Config::rast_jobs)
{
    $csD = "$FIG_Config::rast_jobs/$dataD/CloseStrains";
}

our $dataDF   = "$csD/$dataD";

our $parms = {};

$parms->{genome_types} = &CloseStrains::genome_types($dataDF);
$parms->{dataDF}       = $dataDF;
$parms->{base}         = $base;

if ($request eq "show_otus")
{
    &show_otus($cgi,$csD,$base);
    exit;
}
elsif (($request eq "show_options_for_otu") && $dataD)
{
    $html = &CloseStrains::show_options_for_otu($cgi,$dataD,$dataDF,$base);
    push(@$html,"<a target=_blank href=./wc.cgi?dataD=$dataD&base=$base&request=newick_tree>Newick Tree</a>\n");
}
elsif ($request eq "newick_tree")
{
    if (open(TREE,"<$dataDF/estimated.phylogeny.nwk"))
    {
	my @tree = map { $_ =~ s/,/,\n/g; $_ } <TREE>;
	
	push(@$html,"<pre>",@tree,"</pre>");
	close(TREE);
    }
}
elsif ($request eq "show_signatures")
{
    &CloseStrains::show_signatures($cgi,$parms,$html,$base);
}
elsif ($request eq "show_signatures_reactions")
{
    &CloseStrains::show_signatures_reactions($cgi,$parms,$html,$base);
}
elsif ($request eq "compute_sigs")
{
    &CloseStrains::compute_signatures($cgi,$parms,$html,$base);
}
elsif ($request eq "compute_signatures_reactions")
{
    &CloseStrains::compute_signatures_reactions($cgi,$parms,$html,$base);
}
elsif (($request eq "show_func") && $function)
{
    $function =~ s/^\s+//;
    $function =~ s/\s+$//;
    &CloseStrains::show_func($cgi,$parms,$html,$function,$base);
}
elsif (($request eq "show_family_pegs") && $family)
{
    &CloseStrains::show_family_pegs($cgi,$parms,$html,$family,$base);
}
elsif (($request eq "show_virulence_functions") && (-s "$dataDF/virulence.functions"))
{
    &CloseStrains::show_virulence_functions($cgi,$parms,$html,$base);
}
elsif (($request eq 'show_indexed_funcs') && $keywords)
{
    &CloseStrains::show_indexed_funcs($cgi,$parms,$html,$keywords,$base); 
}
elsif (($request eq "show_ali_or_occurs_tree") && $ali)
{
    &CloseStrains::show_ali($cgi,$parms,$base);  # NOTE: the alignment invokes Gary's alignment viewer,
                                            # which prints the header, so we print everything in show_ali
    exit;
}
elsif (($request eq "show_ali_or_occurs_tree") && $tree)
{
    &CloseStrains::show_occurs_tree($cgi,$parms,$html,$base);
}
elsif (($request eq "show_family_tree") && $family)
{
    &CloseStrains::show_family_tree($cgi,$parms,$html,$family,$base);
}
elsif (($request eq "show_node") && $node)
{
    &CloseStrains::show_node($cgi,$parms,$html,$node,$base);
}
elsif ($request eq "show_otu_tree")
{
    &CloseStrains::show_otu_tree($cgi,$parms,$html,'families',$base);
}
elsif ($request eq "show_reactions")
{
    &CloseStrains::show_reactions($cgi,$parms,$html,$base);
}
elsif (($request eq "show_reaction_on_tree") && $reaction)
{
    &CloseStrains::show_reaction_on_tree($cgi,$reaction,$parms,$html,$base);
}
elsif (($request eq "show_reaction_genome_data") && $reaction && $genome)
{
    &CloseStrains::show_reaction_genome_data($cgi,$reaction,$genome,$dataDF,$parms,$html,$base);
}
elsif ($request eq "show_options_for_reactions")
{
    &CloseStrains::show_options_for_reactions($cgi,$html,$dataD,$dataDF,$base);
#    &CloseStrains::show_signatures_reactions($cgi,$parms,$html,$base);
}
elsif ($request eq "show_otu_tree_reactions")
{
    &CloseStrains::show_otu_tree($cgi,$parms,$html,'reaction',$base);
}
elsif ($request eq "show_adjacency")
{
    &CloseStrains::show_otu_tree($cgi,$parms,$html,'adjacency',$base);
}
elsif ($request eq "show_clusters")
{
    &CloseStrains::show_clusters($cgi,$parms,$html,$base);
}
elsif ($request eq "show_help")
{
    my $type = $cgi->param('type');
    if ($type eq "signatures")
    {
	my $help_sig = &CloseStrains::help_sig();
	push(@$html,"<pre>$help_sig</pre>\n");
    }
    elsif ($type eq "signatures_reactions")
    {
	my $help_sig_reactions = &CloseStrains::help_sig_reactions();
	push(@$html,"<pre>$help_sig_reactions</pre>\n");
    }
    elsif ($type eq "reaction_on_tree")
    {
        my $help_reaction_on_tree = &CloseStrains::help_reaction_on_tree;
        push(@$html,$help_reaction_on_tree);
    }
    elsif ($type eq "compute_signatures")
    {
	my $help_sig_output = &CloseStrains::help_sig_output();
	push(@$html,"<pre>$help_sig_output</pre>\n");
    }
    elsif ($type eq "compute_signatures_reactions")
    {
	my $help_sig_output_reactions = &CloseStrains::help_sig_output_reactions();
	push(@$html,"<pre>$help_sig_output_reactions</pre>\n");
    }
    elsif ($type eq "family")
    {
	my $help_family_output = &CloseStrains::help_family_output();
	push(@$html,"<pre>$help_family_output</pre>\n");
    }
    elsif ($type eq "function")
    {
	my $help_function_output = &CloseStrains::help_function_output();
	push(@$html,"<pre>$help_function_output</pre>\n");
    }
    elsif ($type eq "adjacency")
    {
	my $help_adjacency_output = &CloseStrains::help_adjacency_output();
	push(@$html,"<pre>$help_adjacency_output</pre>\n");
    }
}
else
{
    push(@$html,"<h1>Invalid request</h1>");
}
&HTML::show_page($cgi,$html);
exit;

sub show_otus {
    my($cgi,$datadir,$base) = @_;

    print $cgi->header;
    if (opendir(GENERA,$csD))
    {
	my @genera = grep { $_ !~ /^\./ } readdir(GENERA);
	closedir(GENERA);
	print "<h1>What Changed?</h1>\n";
#	print "<h2><a target='_blank' href=\"http://bioseed.mcs.anl.gov/~overbeek/what_changed.html\">Getting Started: a short Tutorial</a></h2>\n";
	print "<h2>Genera Available</h2>\n";
	foreach my $g (sort @genera)
	{
	    print "<h3><a target='_blank' href='wc.cgi?request=show_options_for_otu&dataD=$g&base=$base'>$g</a>\n";
	}
    }
    else
    {
	print "<h1>The dataD parameter is invalid\n";
    }
}

sub virulence_functions_link {
    my($cgi,$dataDF) = @_;

    if ((-s "$dataDF/virulence.functions") && ($dataDF =~ /([^\/]+)$/))
    {
	my $dataDQ = uri_escape($1);
	return "<a target='_blank' href='wc.cgi?request=show_virulence_functions&base=$base&dataD=$dataDQ'>Some Posssible Functions Associated with Virulence</a>";
    }
    return '';
}
