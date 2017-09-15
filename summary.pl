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

sub random_key {
	my $hash_reference = shift;
	foreach my $key (keys %{ $hash_reference }) {
		return $key;
	}#this subroutine only returns one key of possibly many and then exits
}

sub list_directories {
	my @directories;
	opendir my $dh, '.' or die "Can't opendir on current directory: $!\n";
	while (my $item = readdir $dh) {
		if ($item =~ m/^\.{1,2}$/) {
			next;#skip "." and ".."
		}
		if (-d $item) {
			push @directories, $item;
		}
	}
	closedir $dh;
	return @directories;
}
my %alignment_info;
my %featureCount_info;
foreach my $directory (list_directories()) {
	if ($directory =~ m/Project/) {#I can't access this directory, I don't know why
		next;
	}
	my $key = $directory;
	$key =~ s/_L001$//;
	my $log_file;#this should be the STAR log file
	my $featureCounts_summary;
	chdir "$TOP_DIRECTORY/$directory" or die "Can't chdir to $TOP_DIRECTORY/$directory: $!";
	opendir my $dh, "$TOP_DIRECTORY/$directory" or die "Can't opendir on $TOP_DIRECTORY/$directory: $!";
	while (my $file = readdir $dh) {
		if ($file =~ m/Log\.final\.out$/) {
			$log_file = "$TOP_DIRECTORY/$directory/$file";
		}
		if ($file =~ m/featureCounts\.summary$/) {
			$featureCounts_summary = $file;
		}
	}
	closedir $dh;
	if (!defined $log_file) {
		next;
	}
	open my $log, '<', $log_file or die "Can't read $log_file: $!";
	say "Reading $log_file";
	while (<$log>) {
		if (/^\s+(.+)\s+\|\s+(.+)\s+/) {
			$alignment_info{$1}{$key} =	$2
#			say "\$alignment_info{$1}{$key} =	$2";
		}
	}
	close $log;
#now read featureCounts output
	open my $featureCounts, '<', $featureCounts_summary or die "Can't read $featureCounts_summary: $!";
	while (<$featureCounts>) {
		if (/^Status\s+/) {
			next;
		}
		my @line = split;
		$featureCount_info{$line[0]}{$key} = $line[1];
	}
	close $featureCounts;
	chdir $TOP_DIRECTORY or die "Can't chdir to $TOP_DIRECTORY: $!";
}
my $output_filename = $0;
$output_filename =~ s/\.pl$/.csv/;

my $random_key = random_key(\%alignment_info);
my $hash_key_2 = join ',', sort {$a cmp $b} keys %{ $alignment_info{$random_key} };#make a list of second hash keys

open my $out, '>', $output_filename or die "Can't write $output_filename: $!";
#print all STAR alignment info to the output CSV file
say $out ",$hash_key_2";#this will be the header
foreach my $quality (sort {$a cmp $b} keys %alignment_info) {
	my $hash_key = join ',', sort {$a cmp $b} keys %{ $alignment_info{$quality} };
	if ($hash_key ne $hash_key_2) {
		say "$quality has $hash_key_2, when it should have $hash_key (or maybe vice versa?)";
		die;
	}
	print $out $quality;#this will be the row label
	foreach my $biological_replicate (sort {$a cmp $b} keys %{ $alignment_info{$quality} }) {
		print $out ",$alignment_info{$quality}{$biological_replicate}";
	}
	print $out "\n";#end this row
}
undef %alignment_info;
#print all featureCounts counting to the output CSV
foreach my $quality (sort {$a cmp $b} keys %featureCount_info) {
	my $hash_key = join ',', sort {$a cmp $b} keys %{ $featureCount_info{$quality} };
	if ($hash_key ne $hash_key_2) {
		say "$quality has $hash_key_2, when it should have $hash_key (or maybe vice versa?)";
		die;
	}
	print $out $quality;#this will be the row label
	foreach my $biological_replicate (sort {$a cmp $b} keys %{ $featureCount_info{$quality} }) {
		print $out ",$featureCount_info{$quality}{$biological_replicate}";
	}
	print $out "\n";#end this row
}
close $out;
say "Wrote $output_filename";
