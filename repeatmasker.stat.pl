use strict;

my $f=$ARGV[0];
my %repeat=();
my %Redetail=();
open(IN,"$f")||die "cannot open $f :$!\n";
while(<IN>){
    chomp;
    my @sRows=split/\s+/,$_;
    my $scaf=$sRows[0];
    my $start=$sRows[3];
    my $end=$sRows[4];
    my $strand=$sRows[6];
    #my $repID=$sRows[8];
    my($C,$F,$T)=split/\;/,$sRows[8];
    #$repeat{$repID}+=$end-$start;
    $repeat{$C}+=$end-$start+1;
    if($C ne "Unknown"){
        $Redetail{$C}{$F}+=$end-$start+1;
    }else{
        $Redetail{$C}{$T}+=$end-$start+1;
    }
}
close(IN);

my $out=$f.".stat";
my $out2=$f."detail.stat";
open(OUT,">$out");
foreach my $rep(sort keys %repeat){
    print OUT "$rep\t$repeat{$rep}\n";
}
close(OUT);

open(OUT,">$out2");
foreach my $c(sort keys %Redetail){
    foreach my $f(sort keys %{$Redetail{$c}}){
        print OUT "$c\t$f\t$Redetail{$c}{$f}\n";
    }
}

print "output = $out\n";
print "detail RE = $out2\n";
