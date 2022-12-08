#!/usr/bin/perl -w

#perl AI_index.batch.pl <aa.fasta> <metazoan.dmnd> <NONmetazoan.dmnd>

use strict;
my $dir=$ARGV[0];
my $db_Meta="\/mnt\/nas2\/lihowfun\/bin\/NR_database\/metazoan.dmnd";
my $db_nonMeta="\/mnt\/nas2\/lihowfun\/bin\/NR_database\/NONmetazoan.dmnd";


opendir(Dir,$dir)||die "cannot open $dir:$!\n";
while(my $file=readdir(Dir)){
    if($file=~/(\w+)\.fasta/){
	#### run Diamond blast
	my $command="perl blast_meta.diamond.pl $file $db_Meta $db_nonMeta";
	#system($command);
	
	#### estimate AI index
	my $ID=$1;
	my $meta=$ID."_Meta"; # Diamond blast against Metazoan
	my $nonmeta=$ID."_nonMeta"; # Diamond blast against non-Metazoan
	my $GH=$ID."_out/hmmer.out";
	$command="perl AI_index.pl $meta $nonmeta $GH Orthogroups.txt $file";
	print "$command\n";
	system($command);

		
    }


}


