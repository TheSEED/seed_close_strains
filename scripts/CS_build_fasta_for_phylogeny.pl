use strict;
use Data::Dumper;
use Getopt::Long;
use SeedEnv;
use gjoseqlib;

my $usage = "usage: CS_build_fasta_for_phylogeny -d Data\n";
my $dataD;
my $size = 20000;

my $rc  = GetOptions('d=s' => \$dataD,
		     's=i' => \$size);
if ((! $rc) || (! -d $dataD)) { print STDERR $usage; exit }

my @genomes = map { ($_ =~ /^(\S+)/) ? $1 : () } `cat $dataD/genome.names`;
my $N   = @genomes;
open(FAMS,"cut -f1,3 $dataD/families.good |")
    || die "could not open the families";

mkdir("$dataD/FastaForPhylogeny",0777);
mkdir("$dataD/FastaForPhylogeny/Fasta",0777);
open(INDEX,">","$dataD/FastaForPhylogeny/index") || die "could not open $dataD/FastaForhylogeny/index";

my %to_func = map { ($_ =~ /^\S+\t(\S[^\t]*\S)\t\S+\t(\S+)/) ? ($2 => $1) : () } `cat $dataD/families.all`;
my %to_seq;
foreach my $g (@genomes)
{
    my @seqs    = &gjoseqlib::read_fasta("$dataD/Seqs/$g");
    foreach my $tuple (@seqs)
    {
	my($peg,,undef,$seq) = @$tuple;
	$to_seq{$peg} = $seq;
    }
}

my $n = 1;
my $last = <FAMS>;
while ($last && ($size > 0) && ($last =~ /^(\S[^\t]*\S)\t(\S+)/))
{
    my $func = $1;
    my %pegs;
    my %genomes;
    while ($last && ($last =~ /^(\S[^\t]*\S)\t(\S+)/) && ($1 eq $func))
    {
	my $peg = $2;
	my $g = &SeedUtils::genome_of($peg);
	$genomes{$g} = 1;
	push(@{$pegs{$g}},$peg);
	$last = <FAMS>;
    }

    if (keys(%genomes) == $N)
    {
	my $min = 100000;
	my $max = 0;
	my @filtered_pegs = &filter_pegs(\%pegs,\%to_seq);
	foreach my $peg (@filtered_pegs)
	{
	    my $seq  = $to_seq{$peg};
	    if (length($seq) < $min) { $min = length($seq) }
	    if (length($seq) > $max) { $max = length($seq) }
        }

	if (($max - $min) < 50)
	{
	    print INDEX "$n\t$func\n";
	    open(FASTA,">$dataD/FastaForPhylogeny/Fasta/$n") 
		|| die "cannot open $dataD/FastaForPhylogeny/Fasta/$n";
	    $n++;

	    foreach my $peg (@filtered_pegs)
	    {
		my $g    = &SeedUtils::genome_of($peg);
		my $seq  = $to_seq{$peg};
		print FASTA ">$g $peg\n$seq\n";
	    }
	    $size -= length($to_seq{$filtered_pegs[0]});
	    close(FASTA);
	}
    }
}
close(INDEX);

sub filter_pegs {
    my($pegs,$to_seq) = @_;

    my @keep;
    foreach my $g (keys(%$pegs))
    {
	my $pegs_for_g = $pegs->{$g};
	my @with_ln = sort { $a->[1] <=> $b->[1] } map { [$_,length($to_seq{$_})] } @$pegs_for_g;
	push(@keep,$with_ln[int(@with_ln/2)]);
    }
    return map { $_->[0] } @keep;
}

