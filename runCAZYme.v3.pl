use strict;
# This script help batch running the CAZyme annotation  


my $dir=$ARGV[0];
my %SpeGH=();
my @spearray;

# set the shell script first
# used the conda environment of CAZyme (https://github.com/linnabrown/run_dbcan)
# source activate dbcan

#create a folder and put all the protein fasta file in <specie name>.fasta format 
#batch run_dbcan
opendir(DIR,$dir)||die "cannot open $dir :$!\n";
while(my $file = readdir(DIR)){
    chomp;
    if($file=~/(\S+)\.fasta/){
	my $spe=$1;
	my $output="$spe"."_out";
	push @spearray,$spe;
	my $command = "run_dbcan $spe.fasta protein --db_dir /mnt/nas2/lihowfun/bin/run_dbcan/db.2021 --hmm_cpu 16 --out_dir $output --use_signalP USE_SIGNALP";
	print $command."\n";
	system($command)."\n";
	my $hmmout="$output/hmmer.out";

	#filter the output and count the CAZyme copy numbers of each nematodes 
	my ($length)=&dbcanHmmFilter2($hmmout);
	my %length=%$length;
	foreach my $En(sort keys %length){
	    foreach my $ID(sort keys %{$length{$En}}){
		$SpeGH{$En}{$spe}+=1;		
	    }
	}
    }
}

# print all the copy numbers together
open(OUT,">HmmNum.txt");
print OUT "GH";
for my $i(@spearray){
    print OUT "\t$i";
}
print OUT "\n";
foreach my $GH(sort keys %SpeGH){    
    my $G=$GH;
    substr($G,-4)='';
    print OUT "$G";
    for my $s(@spearray){
	if($SpeGH{$GH}{$s}){
	    print OUT "\t$SpeGH{$GH}{$s}";
	}else{ 
	    print OUT "\t0";
	}
    }
    print OUT "\n";
}


#version 2 (dbcan)
sub dbcanHmmFilter2{
    my $infile=$_[0];
    my %length=();
    my $filtered_file="filtered.".$infile;
    open(OUT,">$filtered_file");
    open(IN,$infile)||die "cannot open $infile :$!\n";
    while(<IN>){
        chomp;
        next if $_=~/\#/;
        my @Rows=split(/\s+/,$_);
        my $GH=$Rows[0];
        my $ID=$Rows[2];
        my $l=$Rows[3];
        my $evalue=$Rows[4];
        my $por=$Rows[9];
        if($por<0.35 || $l<80 || $evalue >1e-15){
            print $GH."\n";
            next;
        }
        print OUT $_."\n";
        if($length{$GH}{$ID}){
            if($length{$GH}{$ID}<$l){                
                $length{$GH}{$ID}=$l;
            }else{
                next;
            }
        }else{
	    $length{$GH}{$ID}=$l;
        }
    }
    close(IN);
    return(\%length);
}
