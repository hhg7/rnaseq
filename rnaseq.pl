#!/usr/bin/env perl

use strict; use warnings; use Cwd; use feature 'say';
my $TOP_DIRECTORY = getcwd();
local $SIG{__WARN__} = sub {#kill the program if there are any warnings
	my $message = shift;
	my $fail_filename = "$TOP_DIRECTORY/$0.FAIL";
	open my $fh, '>', $fail_filename or die "Can't write $fail_filename: $!";
	printf $fh ("$message @ %s\n", getcwd());
	close $fh;
	die "$message\n";
};#http://perlmaven.com/how-to-capture-and-save-warnings-in-perl

say "This script creates STAR commands for a directory of fastq files.  It detects the compression, and tells STAR to align the files with compression in mind.";

sub execute {
	my $command = shift;
	print "Executing Command: $command\n";
	if (system($command) != 0) {
		my $fail_filename = "$TOP_DIRECTORY/$0.fail";
		open my $fh, '>', $fail_filename or die "Can't write $fail_filename: $!";
		print $fh "$command failed.\n";
		close $fh;
		print "$command failed.\n";
		die;
	}
}

sub list_regex_files {
	my @regex;
	while (defined $_[0]) {
		my $regex = shift;
		push @regex, $regex;
	}
	my @files;
	opendir my $dh, '.' or die "Can't opendir on current directory: $!";
	while (my $file = readdir $dh) {
		foreach my $regex (@regex) {
			if ($file =~ m/$regex/) {
				push @files, $file;
			}
		}
	}
	closedir $dh;
	return @files;
}

sub help {
	say "Sample execution: 'perl $0 -g /home/user/directory'";
}

use Getopt::Long 'GetOptions';
my $gtf;
my $debug;
my $genomeDir;

GetOptions ('a=s' => \$gtf,
				'debug|d' => \$debug,
				'genome|g=s' => \$genomeDir
				);#'file' indicates string on command line, '=s' means it must be a string

if (defined $gtf) {
	unless (-e $gtf) {
		say "\n$gtf doesn't exist.\n";
		die;
	}
}

if (!defined $genomeDir) {
	say "A genome directory must be specified for STAR, which can be done with the \"-g\" option.";
	help();
	die;
}

unless (-d $genomeDir) {
	say "$genomeDir isn't a directory.";
	help();
	die;
}

if ($debug) {
	say "Debugging is set to ON.\n";
}

my %files;#if this run is paired-end, I don't want to do runs more than once
foreach my $fastq (list_regex_files('\.fastq', '\.fq')) {
	my $file = "$TOP_DIRECTORY/$fastq";
	my $input = $file;#may later have paired end partner added, which $file won't
	if (defined $files{$file}) {
		next;
	}#if this is a paired end file I've already accounted for, ignore it
	$files{$file} = 1;
	my $prefix;
	if ($file =~ m/$TOP_DIRECTORY\/(.+)\.fastq.bz2/) {
		$prefix = $1;
	} else {
		say "Couldn't find a prefix for $file";
		die;
	}
	if ($file =~ m/_R([12])_001\.fastq/) {#Illumina's paired end file notation
		my $side = $1;
		my $other = $file;
		$prefix =~ s/_R[12]_001//;
		if ($side == 1) {#if this file is the R1 partner
			$other =~ s/_R1_001.fastq/_R2_001.fastq/;
			if (!-e $other) {
				say "The paired end partner for $file ($other) doesn't exist.";
				die;
			}
			$input = "$file $other";
			$files{$other} = 1;
		} else {#i.e. R = 2
			$other =~ s/_R2_002.fastq/_R1_001.fastq/;
			if (!-e $other) {
				say "The paired end partner for $file ($other) doesn't exist.";
				die;
			}
			$input = "$file $other";
			$files{$other} = 1;
		}
	}
	my $compression = '';
	if ($file =~ m/\.bz2$/) {#if the file is bz2 compressed
		$compression = '--readFilesCommand bzcat';#tell STAR the file is bz2 compressed
	} elsif ($file =~ m/\.gz$/) {#if the file is gz compressed
		$compression = '--readFilesCommand zcat';#tell STAR the file is gz compressed
	}
	my $STAR_command = "/usr/bin/time STAR --genomeDir $genomeDir --outFileNamePrefix $prefix $compression --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --readFilesIn $input > $prefix.out 2> $prefix.err";
	mkdir $prefix;
	chdir $prefix or die "Can't chdir to $prefix: $!";
	if ($debug) {#i.e. if $debug isn't defined
		say "Would now run STAR as\n\n$STAR_command";
	} else {
		execute($STAR_command);
	}
	if (defined $gtf) {
		my $bam = $prefix . "Aligned.sortedByCoord.out.bam";
		my $command = "featureCounts -a $gtf -o $prefix.featureCounts $bam";
		if ($debug) {
			say "Would now run featureCounts";
		} else {
			execute($command);
		}
	}
	chdir $TOP_DIRECTORY or die "Can't chdir to $TOP_DIRECTORY: $!";
}
