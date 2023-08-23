
use strict;
use gjoseqlib;
use Carp;
use Data::Dumper;
use Getopt::Long;
use Proc::ParallelLoop;

my $usage = "Usage: pg_build_newick_tree [-t n_threads] -d data_dir \n\n";

my ($help, $dir, $thread);

GetOptions("h|help"     => \$help,
           "d|dir=s"    => \$dir,
           "t|thread=n" => \$thread);

$dir && !$help or die $usage;

$thread ||= 4;

my $sub_dir = "$dir/FastaForPhylogeny";
my $seq_dir = "$sub_dir/Fasta";
my $aln_dir = "$sub_dir/Align"; run("mkdir -p $aln_dir");
my $index_f = "$sub_dir/index";
my $aln_f   = "$dir/concat.aln";
my $tree_f  = "$dir/estimated.phylogeny.nwk";

my @files = map { /(\S+)/ ? $1 : () } `cut -f1 $index_f`;

pareach \@files, sub {
    my ($f) = @_;
    my $cmd = "svr_align_seqs -l -z -muscle < $seq_dir/$f > $aln_dir/$f";
    run($cmd);
}, { Max_Workers => $thread };

my %seqH;
foreach my $f (@files) {
    my @seqs = gjoseqlib::read_fasta("$aln_dir/$f");
    for (@seqs) { $seqH{$_->[0]} .= $_->[2]; }
}

my @merged;
for (sort keys %seqH) {
    push @merged, [ $_, undef, $seqH{$_} ];
}
gjoseqlib::write_fasta($aln_f, @merged);
run("svr_tree -l -c 20 -s < $aln_f > $tree_f");

sub run {
    system($_[0]) == 0 or confess("FAILED: $_[0]");
}
