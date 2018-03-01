#!/usr/bin/perl
# This script searches a .fasta file against a local copy of the swissprot database. Top 10 hits as well as the input (11 total) are printed to hits.txt, their sequences to hits.fasta and a multiple sequence alignment to hits.aln, all stored in a new directory. Finally, hits.aln is opened with ClustalX.
# If no arguments are given, the script does exactly as the assignment requires. The .fasta file to be used may be specified at the command line, or else P11802.fasta is used by default. If the database file is not in the current directory, the script will search for an alternative in another directory.
use strict;
use warnings;
use POSIX 'strftime';
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Perl;
use Bio::Tools::Run::Alignment::Muscle;


# This will be used in naming output folder
my $time = strftime('%Y%m%d_%H_%M_%S', localtime);

# Logfile:
my $LFH;
open ($LFH, ">logfile.txt") or die ("Couldn't create a logfile\n");


# Getting query details, attempting to handle any poor inputs
my $query = 'P11802.fasta';
my $qcheck;

if (@ARGV) {
	$query = $ARGV[0];
	unless (-e $query) {
		$query .= '.fasta';
	}
	unless (-e $query) {
		print ("Query file not found, please give the full name of a .fasta file\n\n");
		print ("Usage: perl bioperl_sjrh2.pl <filename>\n\n");
		print (".fasta files present: \n");
		system ("ls -al *.fasta");
		exit;
	}
}

open (QUERY, "<$query") or die ("Query file not found");

while (<QUERY>) {
	if ($_ =~ /^>sp\|(.+)\|/) {
		$qcheck = $1;
	}
	else {
		&dualprint($LFH, "Query file may be improperly formatted, please use a file from Swissprot if unsatisfactory results are given \n");
	}
	last;
}
close QUERY;

# Locate uniprot_sprot.fasta only if it is not in pwd
my $database = "uniprot_sprot.fasta";
my @locations; 

unless (-e $database) {
	my $uniprots = 	`locate uniprot_sprot.fasta`;
	@locations = split("\n", $uniprots);
	# This loop avoids inappropriate files, e.g uniprot_sprot.fasta.gz	
	foreach (@locations) {
		if ($_ =~ /uniprot_sprot.fasta$/) {
			$database = $_;
			last;
		}
	}
}

if (defined ($locations[1])) {
	&dualprint($LFH, "Multiple files named 'uniprot_sprot.fasta' located. \nUsed for this run: \n\t" . $locations[0] . "\nOthers include:\n");
	foreach (@locations) {
		my $dbloc = "\t" . $_ . "\n";
		&dualprint ($LFH, $dbloc);
	}
}

# Performing BLAST search and outputting results to result.bls
my $blastfac;
my $result;

$blastfac = Bio::Tools::Run::StandAloneBlastPlus->new(
	-db_name => 'swissprot_sjrh2',
	-db_data => $database,
	-create => 1
);

$result = $blastfac->blastp(
	-query => $query,
	-outfile => 'result.bls',
	-method_args => [ '-num_alignments' => 11 ]
);


# Printing BLAST results data to hits.txt
my @IDs;
my $count = 0;

open (HITXT, ">hits.txt") or die ("Couldn't create hits.txt");
print HITXT ("ID\tScore\tE-val\t% Identity\n");

while (my $hit = $result->next_hit) {
	if ($count == 11) {last;}	

	printf HITXT ($hit->accession . "\t" . $hit->bits . "\t" . $hit->significance . "\t%.0f\n", $hit->hsp->percent_identity);
	push (@IDs, $hit->accession);
	$count++;
}

close HITXT;


# Collecting sequences of hits and printing to hits.fasta. $hit->hsp->hit_string would be more efficient, but includes alignment gaps
my $seq_object;

open (HITFA, ">hits.fasta") or die ("Couldn't create hits.fasta");

foreach my $fa (@IDs) {
	$seq_object = get_sequence('swissprot', "$fa");
	print HITFA (">sp|$fa|" . $seq_object->id() . $seq_object->description . "\n" . $seq_object->seq() . "\n");
}

close HITFA;


# Performing multiple sequence alignment and writing results
my $alnfac = Bio::Tools::Run::Alignment::Muscle->new();

my $aln = $alnfac->align('hits.fasta');
my $outaln = Bio::AlignIO->new(
	-format => 'clustalw',
	-file => '>hits.aln'
);
$outaln->write_aln($aln);


# Move results into a separate directory and open alignment in ClustalX
my $dirname = "sjrh2_$qcheck\_$time";
system ("mkdir $dirname");
system ("mv hits.txt $dirname/hits.txt");
system ("mv hits.fasta $dirname/hits.fasta");
system ("mv hits.aln $dirname/hits.aln");
system ("mv result.bls $dirname/results.bls");
system ("mv logfile.txt $dirname/logfile.txt");
system ("clustalx $dirname/hits.aln");


# Tidy up database files. Very specific about filenames so should be safe, but switched off by default.
if (0) {
	system ("rm -f swissprot_sjrh2.phr");
	system ("rm -f swissprot_sjrh2.pin");
	system ("rm -f swissprot_sjrh2.pog");
	system ("rm -f swissprot_sjrh2.psd");
	system ("rm -f swissprot_sjrh2.psi");
	system ("rm -f swissprot_sjrh2.psq");
}


# This subroutine is used in a few places to print the same text to STDOUT and a file
sub dualprint {
	my $fhan = shift;
	my $text = shift;
	my $sout = *STDOUT;

	for my $fh ($fhan, $sout) {
		print $fh ($text);
	}
}
