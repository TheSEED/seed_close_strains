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

# usage: label_all_nodes < tree > relabeled.tree

use tree_utilities;
use Data::Dumper;

($tree = &parse_newick_tree(join("",<STDIN>)))
    || die "could not parse the tree";
&tree_utilities::unlabel_internal_nodes($tree);
&tree_utilities::label_all_nodes($tree);
print &to_newick($tree);
