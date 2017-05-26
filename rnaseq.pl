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

print "This script creates STAR commands for a directory of fastq files.  It detects the compression, and tells STAR to align the files with compression in mind.\n";
print "\n";

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
		say $regex;
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
		print "\n$gtf doesn't exist.\n\n";
		die;
	}
}

if (!defined $genomeDir) {
	say "A genome directory must be specified for STAR, which can be done with the \"-g\" option.";
	die;
}

unless (-d $genomeDir) {
	say "$genomeDir isn't a directory.";
	die;
}

if ($debug) {
	say "Debugging is set to ON.\n";
}

foreach my $file (list_regex_files('\.fastq', '\.fq')) {
	print "$file\n";
	my $prefix;
	if ($file =~ m/(\d+)\.fastq.bz2/) {
		$prefix = $1;
	} else {
		print "Couldn't find a prefix for $file\n";
		die;
	}
	my $compression = '';
	if ($file =~ m/\.bz2$/) {#if the file is bz2 compressed
		$compression = '--readFilesCommand bzcat';#tell star the file is bz2 compressed
	} elsif ($file =~ m/\.gz$/) {#if the file is gz compressed
		$compression = '--readFilesCommand zcat';#tell star the file is gz compressed
	}
	my $STAR_command = "STAR --genomeDir $genomeDir --outFileNamePrefix $prefix $compression --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40 --readFilesIn $TOP_DIRECTORY/$file > $prefix.out 2> prefix.err";
	mkdir $prefix;
	chdir $prefix or die "Can't chdir to $prefix: $!";
	if (!$debug) {#i.e. if $debug isn't defined
		print "Would now run STAR\n";
#		execute($STAR_command);
	}
#	chdir '..' or die "Can't chdir to original directory: $!";
	if (defined $gtf) {
		my $bam = $prefix . "Aligned.sortedByCoord.out.bam";
		my $command = "featureCounts -a $gtf -o $prefix.featureCounts $bam";
		if (!$debug) {
			print "Would now run featureCounts\n";
#			execute($command);
		}
	}
}

