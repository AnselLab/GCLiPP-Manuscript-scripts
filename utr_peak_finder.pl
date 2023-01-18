#Creator: Adam Litterman
#Go through a list of bed regions
#and find peaks of a specified width in each one

#Works by sliding a normal distribution with a height of depth(X) along 
#the entire UTR and finding the correlation between that and the observed depth surrounding the 
#location. Then goes through and finds local maxima in the correlation

#This script updates peak_finder_4 by working with different scoring algorithms to report peaks 
#based on several criteria (not just depth of protein depth file)
#and also having an option to report only the n best scoring peaks

#USAGE EXAMPLES:
#perl utr_peak_finder.pl -debugUTR -gene=Il2ra
#perl utr_peak_finder.pl -bed=refseq_mm10_coding_3putrs_extended_092815.bed -printall > test_peaks_092815.bed
#perl utr_peak_finder.pl -mode=inverse_conservation -bed=refseq_mm10_coding_3putrs_extended_092815.bed -printall
#perl utr_peak_finder.pl -mode=inverse_conservation -bed=refseq_mm10_coding_3putrs_extended_092815.bed -n=5000
#perl utr_peak_finder.pl -mode=scaled -bed=refseq_mm10_coding_3putrs_extended_092815.bed -printall
#perl utr_peak_finder.pl -mode=inverse_scaled -printall


use List::Util qw(sum min max);

use Data::Dumper qw(Dumper);


foreach $parameter( @ARGV )
{
	if( substr($parameter,0,1) eq '-')
	{
		($pname, $pval) = split( '=',substr($parameter,1)) ;
		if ($pval eq '') { $pval = 1; }
		$PARAMS{$pname} = $pval;
	}
		
}


$utrfile = "/Users/DarthRNA/Documents/Adam/genomes/refseq_mm10_custom_3putr_comprehensive_092115.bed";


$utrfile = "/Users/DarthRNA/Documents/Adam/genomes/custom_short.bed" if $PARAMS{short};
$utrfile = $PARAMS{bed} if $PARAMS{bed};


#$utrfile = "/Users/DarthRNA/Documents/Adam/genomes/refseq_mm10_custom_3putr_comprehensive_092115.bed";
$symbolfile = "/Users/DarthRNA/Documents/Adam/genomes/refseq_symbols.txt";


#CONSTANT parameters for normal mode peak calling
$MINIMUM_SMALLEST_ACCEPTABLE_PEAK = 12; #constant that is the minimunum number of reads that a peak must be in a low expressed gene 

$WINDOW_WIDTH = 71;			   #width of the peaks; should be odd
$DIFF_THRESHOLD = 0.10;        #minimum height of peaks relative to maximum UTR depth


#assign the parameters from runtime command line input if the command line input exists

$MINIMUM_SMALLEST_ACCEPTABLE_PEAK = $PARAMS{minsap} if $PARAMS{minsap}; 

$WINDOW_WIDTH = $PARAMS{ww} if $PARAMS{ww};
$DIFF_THRESHOLD = $PARAMS{dt} if $PARAMS{dt};

$DONTWRITE = 1 if $PARAMS{dontwrite};
$DISCARD_BLOCKS = 0 if $PARAMS{keepblocks};

$HALF = int($WINDOW_WIDTH / 2);
		 

open(SYMBOLFILE, $symbolfile);
while(<SYMBOLFILE>)
{
	chomp();
	($refseqname, $symbol) = split('	');
	$SYMBOLS{$refseqname} = $symbol;
}
close SYMBOLFILE;	


#main loop
my @output_stack = ();
my @score_stack = ();
my @score_stack = ();
my $global_peak_counter = 1;


@utrs = ();

open(UTR, $utrfile);
while(my $utrline = <UTR>)
{
	chomp($utrline);
	my ($chr, $start, $end, $name, $score, $strand) = split('	', $utrline);

	#next if $score < 1 && !$PARAMS{dontskip};

	my $refseq = undef;
	$refseq =  "NM_$1" if $name =~ /NM_(\d+)/;

	$symbol = $SYMBOLS{$refseq};
	$symbol = $name if $refseq eq '' || $symbol eq '';

	my $locationName = $chr . ':' . $start . '-' . $end;

	my $which_bam;
	my $max_utr_depth = 0;
	
	my $prop_ref; 
	my @prop = (); 
	
	($prop_ref, $max_utr_depth) = getProportionalCoverageAndMax($strand, $chr, $start, $end);

	if ($PARAMS{debugUTR} && $PARAMS{gene} eq $symbol)
	{
		print "Index\tNormal\tInverse\tScaled_Normal\tScaled_Inverse\tScaled_Cons_Normal\tScaled_Cons_Inverse\n";
		dumpArrays($prop_ref, $inverse_ref, $scaled_prop_ref, $scaled_inverse_ref, $scaled_cons_prop_ref, $scaled_cons_inverse_ref);
	}

	next unless defined $prop_ref;

	@prop = @{$prop_ref};
	
	if($max_utr_depth > $MINIMUM_SMALLEST_ACCEPTABLE_PEAK)
	{
		
		$utrsize = scalar(@prop);

		my @minima;
  		my @maxima;
  		my @peak_locations;
  		my @heights;
  		
  		my @test_scores;
  		
  	    my $prev_cmp = 0;

		for($i = 0; $i < $utrsize; $i++)
		{	
			if($prop[$i])
			{
				my @peak = ();
				for $q ($i - $HALF .. $i + $HALF)
				{
					if($prop[$q]) {push(@peak, $prop[$q]);}
					else{ push(@peak,0); } 
				}
				my @idealized_peak = getNormal($WINDOW_WIDTH, $prop[$i]);
			
				$score = pearson(\@peak, \@idealized_peak);
				$test_scores[$i] = $score;
				$p = $prop[$i];

			}	
    	}

		for($i = 0; $i < $utrsize; $i++)
		{	
			#see if this entry is bigger or smaller than the last entry
		 	my $cmp = $test_scores[$i] <=> $test_scores[$i+1];
    		 	
    		if ($cmp != $prev_cmp)  #if this comparison is not equal to the previous one...
    		{
    			my $left = 0, $right = 0;
    			$left = $i-$HALF;
    			$left = 0 if $left < 0;
    			$right = $i+$HALF;
    			$right = $#inv if $right > $#inv;

    			my $height = $prop[$i] * $max_utr_depth; 
    				
				my $test = $prop[$i] > $DIFF_THRESHOLD && $prop[$i] && $height > $MINIMUM_SMALLEST_ACCEPTABLE_PEAK;
    			
    			if($cmp > 0 && $test)  #and this entry is bigger than the next one...
    			{
    				#then we've found a local maximum so add it to the list of peaks unless it's very close to the previous
    				#peak in which case we only keep the better of the two
    				$last_maximum = $maxima[-1];
    				$last_location = $peak_locations[-1];
    				$maximum = $test_scores[$i];
    				$peak_location = $i; 
 					
    				#if the peaks are not too close together add the new one
    				if(abs($peak_location - $last_location) > $HALF)
    				{
    					push(@maxima, $maximum);
     					push(@peak_locations, $peak_location);
     					push(@heights, $height);
     				} #if they are, but the new peak is better, then ditch the last one and add the new one
     				elsif( $maximum > $last_maximum)
     				{
     					pop(@maxima);
     					pop(@peak_locations);
     					pop(@heights);
     					push(@maxima, $maximum);
     					push(@peak_locations, $peak_location);
     					push(@heights, $height);
					}
					else #if they are too close but the last one was better do nothing, skip this one
					{}
     						
     			}
     			
     			# when this and next elements are ==, defer checking for
      			# minima/maxima till next loop iteration
      			$prev_cmp = $cmp if $cmp;
    		}
			
    	}
				
		for($m = 0; $m < scalar(@maxima); $m++)
		{	
			my $name =~ s/\s+//g;
			my $peakStart = $peak_locations[$m] - $HALF + $start;
			my $peakEnd = $peak_locations[$m] + $HALF + $start;
			my $peakScore = $heights[$m];
			my $peakLocation = $chr . ':' . $peakStart . '-' . $peakEnd;
	
			my $prefix = '';
			$prefix = $PARAMS{prefix} . '_' if $PARAMS{prefix};
			$peakName = $prefix . $symbol . "_" . $peakLocation . "_" . $global_peak_counter++;

			#$peakName = $PARAMS{prefix} . "_" .  $symbol . "_" . $peakLocation . "_" . sprintf("%x",$global_peak_counter++);
			#$peakName = $symbol . "_" . sprintf("%x",$global_peak_counter++); 

			#print "$chr\t$peakStart\t$peakEnd\t$peakName\t0\t$strand\n" if $DONTWRITE;
			#print "$chr\t$peakStart\t$peakEnd\t$peakName\t$peakHeight\t$strand\n" unless $DONTWRITE;
			
			my $lastline = $output_stack[-1];
			chomp($lastline);
			my ($lastChr, $lastPeakStart, $lastPeakEnd, undef, $lastScore, $lastStrand ) = split("\t",$lastline);

			#print "$lastline $lastChr eq $chr && $lastStrand eq $strand && abs($peakStart - $lastPeakStart) < $HALF\n";

			#if this peak is too close to the last peak
			if( $lastChr eq $chr && $lastStrand eq $strand && abs($peakStart - $lastPeakStart) < $HALF)
			{
				#then pop the last one if this one is better
				if($peakScore > $lastScore)
				{
					pop(@output_stack);
					pop(@score_stack);
					push(@output_stack,"$chr\t$peakStart\t$peakEnd\t$peakName\t$peakScore\t$strand\n");
					push(@score_stack, $peakScore);
				}
				#otherwise do nothing

			}#if not then just add it to the stack
			else
			{
				push(@output_stack,"$chr\t$peakStart\t$peakEnd\t$peakName\t$peakScore\t$strand\n");
				push(@score_stack, $peakScore);
			}
				
		} 	
	}
}

close(UTR);

#print the peaks (if any) for the last UTR
if($PARAMS{printall})
{
	print foreach @output_stack;
}

if($PARAMS{n})
{
	my @sorted_indices = sort { $score_stack[$b] <=> $score_stack[$a] } 0 .. $#score_stack;
	for$c(0..$PARAMS{n})
	{
		print $output_stack[$sorted_indices[$c]];
	}
}


sub getProportionalCoverageAndMax
{
	my ($strand, $chr, $start, $end) = @_;


	#peak calling is done with same-stranded reads so need to sort out strandedness before 
	#generate these files with commands like this
	#samtools view  -b -f 16 jurkat_all.bam > jurkat_all_CRICK.bam 
	#samtools view -b -F 16 jurkat_all.bam > jurkat_all_WATSON.bam
	my $CRICK_bam_file = "all_CRICK.bam";
	my $WATSON_bam_file = "all_WATSON.bam";

#	if($PARAMS{bam})
#	{
#		$CRICK_bam_file = "../alignments/uj_bam/" . $PARAMS{bam} . "_CRICK.bam";
#		$WATSON_bam_file = "../alignments/uj_bam/" . $PARAMS{bam} . "_WATSON.bam";		
#	}

	my $bam;
	$bam = $WATSON_bam_file if $strand eq '+';
	$bam = $CRICK_bam_file if $strand eq '-';

	my @prop = ();

	my $width = abs($end-$start);

	my @prop = (0) x $width;
	my $max = 0;

	#print "samtools depth -r $chr:$start-$end  $bam |\n";
	(open(DEPTH, "samtools depth -r $chr:$start-$end  $bam |")) || die "Couldn't call samtools";	
	while (<DEPTH>) 
	{
		chomp;
		my (undef, $position, $depth) = split ("\t");
		if ($depth>$max){$max = $depth;}

		my $c=$position - $start;
		$prop[$c] = $depth;
	}
	close DEPTH;

	return(undef)x2 if $max == 0;
	$_ /= $max for(@prop);
	return(\@prop, $max);
}	

	
sub getNormal
{
	#getNormal(51, 0.8) 

	#return numerical approximation of a normal distribution from -1 to 1 with sigma^2 = 0.2 as an array of width 51 and height 0.8
	  
	#first number should be odd -- will return an odd number of entries even if first number is even
	  
	my $w = $_[0];
	my $h = $_[1];
	
	my @toreturn = ();

	my $half = int($w/2);
	  
	my $mu = 0;
	my $sigma = 0.4472135955;
	my $pi = 3.14159265359;
	my $e = 2.718281828459;

	
	my $sd = $_[0] / (2*$_[1]);
	
	for $i(-$half .. $half)
	{
		
		$x = $i/$half;
		$b = $e ** ( -2.5 * ( ($x * $x) ) ) ;	 #2.5 is hardcoded b/c sigma^2 = 0.2 
		$q = $i + $half;
        $r = $h * $b;
		$toreturn[$q] = $r;
	}
	
	return @toreturn;
}
	

sub pearson
{

	my @x = @{$_[0]};
	my @y = @{$_[1]};
	
	
	$x_m = A(@x);
	$y_m = A(@y);
	
	$sigma_x = sum(@x);
	$sigma_y = sum(@y);

	$sigma_xy = 0;
	$sigma_sx = 0;
	$sigma_sy = 0; 
	
	$n = scalar @x;
	$len = scalar @x - 1;
	for $i (0 .. $len)
	{
		$sigma_xy += $x[$i] * $y[$i];
		$sigma_sx += ($x[$i] - $x_m) * ($x[$i] - $x_m); 
		$sigma_sy += ($y[$i] - $y_m) * ($y[$i] - $y_m); 
		
	}
	

	$sx = sqrt( (1 / $len) * $sigma_sx );
	$sy = sqrt( (1 / $len) * $sigma_sy );
	
	$top = ( $sigma_xy - ( $n * $x_m * $y_m ) ) ;
	$bottom = ($len * $sx * $sy );
	
	if($bottom != 0) {	return $top / $bottom; }
	else{ return 1; } 
	
} 

sub dumpArrays
{
	my @first = @{$_[0]};
	my @refs = @_;
	my $len = scalar(@first);
	my @line = ();
	for $x(0..$len-1)
	{
		print $x."\t".$_[0]->[$x]."\t".$_[1]->[$x]."\t".$_[2]->[$x]."\t".$_[3]->[$x]."\t".$_[4]->[$x]."\t".$_[5]->[$x]."\n";
		#push(@line, $x);
		#for $p(0..scalar(@_)-1){push(@line,$refs[$p]->[$x])}
		#print $x . "\t" . scalar(@_) . "\t" . join(@line,"\t") . "\n";
	}
}

sub A
{
        my $a = 0;
        $a += $_ for @_;
        return $a / @_;
}

	
