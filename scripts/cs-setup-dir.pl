use strict;
use Data::Dumper;
use CloseStrains;
use File::Path 'make_path';
use File::Basename;
use Getopt::Long::Descriptive;
use POSIX;
use File::Slurp;
use Job48;
use FIGM;

my($opt, $usage) = describe_options("%c %o rast-job-dir extra-genomes [...]",
				    ["help|h" => "Show this help message"],
				   );
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV < 1;

my $dir = shift;
my @extra = @ARGV;

-d $dir or die "Directory $dir does not exist\n";

my $cs_dir = "$dir/CloseStrains";
make_path($cs_dir);

my $job_id = basename($dir);

my $path = "$cs_dir/$job_id";

if (-d $path)
{
    my $now = strftime("%Y-%m-%d-%H-%M-%S", localtime);
    my $bak = "$path.bak.$now";
    rename($path, $bak);
    write_file("$bak/HIDE", "Hidden due to creation of new set\n");
}
make_path($path);

my $job = Job48->new($dir);
my $fig = FIGM->new(undef, $job->orgdir);

my $gobj = $fig->genome_id_to_genome_object($job->genome_id);

my @rast_genomes;
my @refs;

push(@rast_genomes, [$job_id, $gobj]);

for my $e (@extra)
{
    # if ($e =~ /^\d+$/)
    # {
    # 	my $ejob = $self->app->data_handle('RAST')->Job->init({ id => $e });
    # 	if ($ejob)
    # 	{
    # 	    my $gobj = $self->job_to_genome_object($ejob);
    # 	    if ($gobj)
    # 	    {
    # 		push(@rast_genomes, [$ejob->id, $gobj]);
    # 	    }
    # 	    else
    # 	    {
    # 		warn "Failed to create gobj for job $e\n";
    # 	    }
    # 	}
    # 	else
    # 	{
    # 	    warn "Cannot find RAST job for job id $e\n";
    # 	}
    # }
    if ($e =~ /^(core|p3)\|\d+\.\d+$/)
    {
	#
	# This is a PATRIC or coreseed genome. Save it in the refs list for later expansion.
	#
	push(@refs, $e);
    }
    elsif ($e =~ /^\d+\.\d+$/)
    {
	#
	# We check first to see if this is a valid genome in pubseed
	# (since all PubSEED genomes went thru rast, we prefer the pubseed
	# version to the RAST version).
	#
	my $gdir = "/vol/public-pseed/FIGdisk/FIG/Data/Organisms/$e";
	if (-d $gdir && ! -f "$gdir/DELETED")
	{
	    push(@refs, $e);
	}
	# else
	# {
	#     my $ejob = $self->app->data_handle('RAST')->Job->init({ genome_id => $e });
	#     if ($ejob)
	#     {
	# 	my $gobj = $self->job_to_genome_object($ejob);
	# 	if ($gobj)
	# 	{
	# 	    push(@rast_genomes, [$ejob->id, $gobj]);
	# 	}
	# 	else
	# 	{
	# 	    warn "Failed to create gobj for RAST genome $e\n";
	# 	}
	#     }
	#     else
	#     {
	# 	push(@refs, $e);
	#     }
	# }
    }
    else
    {
	warn "Unparsable job specfifier '$e'\n";
    }
}

CloseStrains::create_set_from_rast($path, \@refs, \@rast_genomes);
CloseStrains::get_genome_name($path);
