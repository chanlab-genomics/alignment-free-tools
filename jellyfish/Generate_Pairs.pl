use strict;
my @name;
open (TAXA, "files2run.txt");
while (<TAXA>) {
	chomp $_;
	my @array = split ("\t", $_);
	push (@name, $array[0]);
	}
close (TAXA);

my %check; my $k = 1;
for (my $i = 0; $i < scalar(@name); $i++) {
	for (my $j = $k; $j < scalar(@name); $j++) {
		print $name[$i]."\t".$name[$j]."\n";
		}
	$k++;
	}
