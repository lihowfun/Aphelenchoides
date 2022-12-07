use strict;
my $inrepgff=$ARGV[0];
my $outfile="categorize.".$inrepgff;
open(OUT,">$outfile");
my %repeat=();
open(IN,$inrepgff);
while(<IN>){
    chomp $_;
    next if $_=~/query/;
    next if $_=~/sequence/;
   # next if $_=~/RepeatMasker/;
    #next if $_=~/^\s+/;
    $_=~s/^\s+//;
    my @sRows=split/\s+/,$_;
    next unless $sRows[10] ne '';
    my($C,$F)=split/\//,$sRows[10];
    $C=~s/\?//g;
    if($sRows[8] eq '+'){
        if($F){
            print OUT "$sRows[4]\tRepeatMasker\tsimilarity\t$sRows[5]\t$sRows[6]\t$sRows[1]\t+\t.\t$C;$F;$sRows[9]\n";
        }else{
            print OUT "$sRows[4]\tRepeatMasker\tsimilarity\t$sRows[5]\t$sRows[6]\t$sRows[1]\t+\t.\t$C;NA;$sRows[9]\n";
        }
    }else{
        if($F){
            print OUT "$sRows[4]\tRepeatMasker\tsimilarity\t$sRows[5]\t$sRows[6]\t$sRows[1]\t-\t.\t$C;$F;$sRows[9]\n";

        }else{
            print OUT "$sRows[4]\tRepeatMasker\tsimilarity\t$sRows[5]\t$sRows[6]\t$sRows[1]\t-\t.\t$C;NA;$sRows[9]\n";
        }
    }
}

print "my output file = $outfile\n";
