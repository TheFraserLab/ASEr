#!/usr/bin/perl


############################################
# UPDATE HISTORY AND LIST OF THINGS TO FIX #
############################################

#12.09.26 - Initial work on script.


###########################
# SUBROUTINES AND MODULES #
###########################

#Use the module that allows for long command line option parsing. 
use Getopt::Long;
use Term::ANSIColor;

sub fasta2hash	{
	my($file) = @_; 
	my @FASTA;
	my $line;
	my @line2;
	my @line3;
	my %fastahash;
	my $curhead;
	my $seq = "";
	my $tmp;

	#Open the FASTA file and store it in an array.

	open (LIST, "$file");
	@FASTA = <LIST>;
	close LIST;

	#Now go through the FASTA an store '>' lines as KEYS and sequence as VALUES.

	foreach $line (@FASTA)	{
		chomp($line);	
	
		if(($line =~ />/) && ($seq eq ""))	{	#Here's what we do with the first header.
			@line2 = split(/>/, $line);
			@line3 = split(/ /, $line2[1]);
			$curhead = $line3[0];
		}

		if($line !~ />/)	{	#Here's what we do with sequence lines.
			@line2 = split(/>/, $line); 
			$seq .= $line;
		}
		
		if(($line =~ />/) && ($seq ne ""))	{	#Here's what we do with the subsequent headers.
			$fastahash{$curhead} = $seq;
			$seq = "";
			@line2 = split(/>/, $line);
			@line3 = split(/ /, $line2[1]);
			$curhead = $line3[0];
		}
	}
	
	$fastahash{$curhead} = $seq;	#The final FASTA seq will be put in here.
	
	return %fastahash;
}

if(($ARGV[0] eq '-h') || ($ARGV[0] eq '--help'))	{
	print colored['bright_red'], '
	This script takes a list of SNPs in BED format as well as a genome in FASTA format and 
	outputs a a new FASTA file containing the genome with all SNP positions masked as \'N\'s.

	USAGE: MaskReferencefromBED.pl <SNP BED FILE> <GENOME FASTA FILE> <MASKED OUTPUT FASTA>
	
		A list of SNPs in BED format must be supplied as follows:
   
		CHR \t 0-POSITION \t 1-POSITION \t REF|ALT
   
		e.g.
   
		chr02	1242	1243	A|G

	--help or -h
		Print this text.
		
';
	exit;
}


######################################
# DEFAULT AND COMMAND LINE VARIABLES #
######################################

##########
# SCRIPT #
##########

#Read in BED file
open(BED, "$ARGV[0]");
while(<BED>)	{
	@line = split(/\t/, $_);
	$pos = "$line[0]\_$line[1]";
	$bedhash{$line[0]} .= "$line[1]\t";
	$chrs{$line[0]} = 1;
}
close BED;

#Read in the genome
%genomehash = fasta2hash($ARGV[1]);

@allchrs = keys(%genomehash);

#How many chromosomes have variants?
@keys = keys(%chrs);
$varchrs = scalar(@keys);
$totchrs = scalar(@allchrs);
print "$varchrs OUT OF $totchrs CHROMOSOMES/CONTIGS HAVE AT LEAST ONE VARIANT\n";

open(OUT, ">$ARGV[2]");
$maskedsites = 0;
foreach $chr (@allchrs)	{
	@seq = split('', $genomehash{$chr});
	
	@mask = split(/\t/, $bedhash{$chr});
	
	for $m (@mask)	{
		$seq[$m] = "N";
		++$maskedsites; 
	}
	$outseq = join('', @seq);
	print OUT ">$chr\n$outseq\n\n";
}
close OUT;

print "$maskedsites SITES MASKED\n\n";
