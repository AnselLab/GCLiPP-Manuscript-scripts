Creator: Adam Litterman
use List::Util 'max';
use List::Util 'sum';
use Statistics::Basic qw(:all);
use Data::Dumper qw(Dumper);
use Statistics::Distributions;


foreach $parameter( @ARGV )
{
	if( substr($parameter,0,1) eq '-')
	{
		($pname, $pval) = split( '=',substr($parameter,1)) ;
		if ($pval eq '') { $pval = 1; }
		$PARAMS{$pname} = $pval;
	}
		
}

$N = 2;

#number of species for which there must be sequences to align
$N = $PARAMS{n} if $PARAMS{n};

#Print 3' UTRs for which there is one and only one entry corresponding to a given symbol (to avoid difficulty of aligning alternatively spliced 3' UTRs)
my %symbols = {};
my %bed_entries = {};

#my $fasta_dir = '/Users/DarthRNA/Documents/Frisco/phylogeny_aug16/fasta';
#my $gene_dir = '/Users/DarthRNA/Documents/Frisco/phylogeny_aug16/genes';
#my $alignment_dir = '/Users/DarthRNA/Documents/Frisco/phylogeny_aug16/alignments';

#my $fasta_dir = '/Users/DarthRNA/Documents/Frisco/phylogeny_jan18/fasta';
#my $gene_dir = '/Users/DarthRNA/Documents/Frisco/phylogeny_jan18/genes';
#my $alignment_dir = '/Users/DarthRNA/Documents/Frisco/phylogeny_jan18/alignments';

my $fasta_dir = '/Users/DarthRNA/Documents/Frisco/compare_jan18/peaks/peak_alignment/fasta';
my $gene_dir = '/Users/DarthRNA/Documents/Frisco/compare_jan18/peaks/peak_alignment/regions';
my $alignment_dir = '/Users/DarthRNA/Documents/Frisco/compare_jan18/peaks/peak_alignment/alignments';



my %sequences = {};

opendir(my $fh, $fasta_dir) || die "Can't open $fasta_dir: $!";
while (readdir $fh) 
{
	next if $_ eq '.' || $_ eq '..';
	my $fasta_file = $fasta_dir . '/' . $_;
	my ($organism, undef) = split('_'); 
	open F, $fasta_file;
	while(my $symbol_line = <F>)
	{

		my $sequence = <F>;
		chomp $symbol_line;
		chomp $sequence;
		my $symbol = substr($symbol_line,1);

		$sequences{$symbol}{$organism} = $sequence;
		$sequences{$symbol}{count}++; 
	}
	close F;

}
closedir $fh;


#print Dumper \%sequences;
@counts = (0) x 10;


foreach $symbol (keys %sequences)
{
	$counts[$sequences{$symbol}{count}]++;

	next unless $sequences{$symbol}{count} >= $N;
	my $gene_file = $gene_dir . '/' . $symbol . '.fa';
	open G, ">$gene_file";
	foreach $species ( sort{ $a cmp $b } keys %{$sequences{$symbol}} )
	{
		next if $species eq 'count';
		print G ">$species\n$sequences{$symbol}{$species}\n";
	}
	close G;

}


opendir(my $gh, $gene_dir) || die "Can't open $gene_dir: $!";
while (readdir $gh) 
{
	next if $_ eq '.' || $_ eq '..';
	my $gene = substr($_,0,-3);
	my $gene_file = $gene_dir . '/' . $_;
	my $alignment_file = $alignment_dir . '/' . $gene . '_aligned.fa';
	print "clustalo -i $gene_file -o $alignment_file -v\n";
	open (C, "clustalo -i $gene_file -o $alignment_file|") || die;
	close C;
}
closedir $gh;
#print Dumper \@counts;
