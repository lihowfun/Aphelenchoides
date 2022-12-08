use strict;

my $aa=$ARGV[0];

#data base location, should be modified accordingly
my $db_Meta="\/mnt\/nas2\/lihowfun\/bin\/NR_database\/metazoan.dmnd";
my $db_nonMeta="\/mnt\/nas2\/lihowfun\/bin\/NR_database\/NONmetazoan.dmnd";

my $spe_ID;
if($aa=~/(\w+)\.fasta/){
    $spe_ID=$1;
}

my $Meta_out="$spe_ID"."_Meta";
my $nonMeta_out="$spe_ID"."_nonMeta";
my $NR_out="$spe_ID"."_NR";

if(-e $Meta_out){
}else{
    &run_diamond($aa,$db_Meta,"$Meta_out",1);
}

if(-e $nonMeta_out){
}else{
    &run_diamond($aa,$db_nonMeta,"$nonMeta_out",1);
}


#diamond latest version is v2.0.x
sub run_diamond{
    my ($fasta,$db,$out,$target_Num)=@_;
    #my $db=$_[1];
    print $db."!!!!\n";
    #my $out=$_[2];
        my $command=" diamond blastp \\
        --query $fasta \\
        --max-target-seqs $target_Num \\
        --sensitive \\
        --threads 32 \\
        --db $db \\
        --outfmt 6 \\
        --out $out  \\
        --evalue 1e-3 ";
    print "$command\n";
    system($command);
}

