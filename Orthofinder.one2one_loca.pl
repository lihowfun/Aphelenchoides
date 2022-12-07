use strict;

my $Ingroup=$ARGV[0];
my $Insequence=$ARGV[1];
my $Inspe=$ARGV[2];
my $spe1gff=$ARGV[3];
my $spe2gff=$ARGV[4];
#my $proteinID=$ARGV[5];
my $target_family=$ARGV[6];


    
my @sRow0='';
my @species='';

my %protein=();
my %geneID=();

my @sRows;
my $gene;

my %chr=();
my %start=();
my %end=();
my %strand=();
my $spe1=$spe1gff;
substr($spe1,-4)='';
my $spe2=$spe2gff;
substr($spe2,-4)='';


print "Usage = perl GetOne2oneLoc.v1.1.pl Orthogroups.txt SequenceIDs.txt SpeciesIDs.txt Chicken.gff Mandarin.gff proteinID_list <OG0000003>\n";

#JGI A.oligospora
open(IN,$Inspe)||die "cannot open $Inspe: $!";
while(<IN>){
    chomp $_;
    @sRow0=split/\s+/,$_;
    substr($sRow0[1],-6)='';
    push @species,$sRow0[1];
    
}
close(IN);

foreach my $sp(@species){
    print $sp."\n";
}
print "\n";

open(IN,$spe1gff)||die "cannot open $spe1gff $! ";
while(<IN>){
    next unless ($_ =~ /mRNA/ && $_ =~ /Parent/);
    next if ($_ =~ /CDS/);
    chomp $_;
    if($_=~/ID=(\w+):/){
	my $geneID=$1;
	@sRows=split /\s+/,$_;
	$chr{$geneID}=$sRows[0];
	# print $chr{$geneID}."\n";
	$start{$geneID}=$sRows[3];
	$end{$geneID}=$sRows[4];
	$strand{$geneID}=$sRows[6];
    }
}
close(IN);

#Maker
@sRows='';
open(IN,$spe2gff)||die "cannot open $spe2gff $! ";
while(<IN>){
    next unless ($_ =~ /mRNA/ && $_ =~ /Parent/);
    next if ($_ =~ /CDS/);
    chomp $_;
    #print $_."\n";
    if ($_=~/ID=(\w+):/){
	my $geneID=$1;
	#print $geneID."\n";
	@sRows=split /\s+/,$_;
	$chr{$geneID}=$sRows[0];
	$start{$geneID}=$sRows[3];
	$end{$geneID}=$sRows[4];
	$strand{$geneID}=$sRows[6];
    }
    
}


@sRows='';
open(IN,$Insequence)||die "cannot open $Insequence: $! ";
while(<IN>){
    chomp $_;
    @sRows=split/\s+/,$_;
    $gene=$sRows[0];
    substr($gene,-1)='';
    $geneID{$gene}=$sRows[1];
}
close(IN);


my @sRows1;
my %count;
my %spegen=();
my %spegen2=();
my @targetA='';
my @targetB='';
open(IN,$Ingroup)||die "cannot open $Ingroup: $! ";
while(<IN>){
    chomp $_;
    @sRows1=split/\s+/,$_;
    my $famID=shift @sRows1;
    substr($famID,-1)='';
    my $v='';
    foreach $v (@sRows1){
	my ($spe,$gene)=split/\|/,$v;

	
	#if($spe eq $spe1){
	 #   substr($gene,-2)='';
	    # print $gene."\n";
	#}
	$count{$famID}{$spe}+=1;
	#$spegen{$gene}=$spe;
	################### need to be modify if we want to parsing 1 genefamily contain mutiple genes in same specie 
	$spegen{$famID}{$spe}=$gene;
	$spegen2{$famID}{$gene}=$spe;
	#print $spegen{$famID}{$spe}."\n";
    }
    
    
    
    #print $famID;
}
close(IN);


print "$spe1\t$spe2\n";

if($target_family){
    my $spe1_target=$spe1."_".$target_family.".txt";
    my $spe2_target=$spe2."_".$target_family.".txt";
    open(tarA,">$spe1_target");
    open(tarB,">$spe2_target");

    foreach my $key1(sort keys %spegen2){
	#print $key1."\n";
	if($key1 eq $target_family){
	    foreach my $gene(sort keys %{$spegen2{$key1}}){
		#print $gene."\n";
		if($spegen2{$key1}{$gene} eq $spe1){
		    #print $gene."\n";
		    if($chr{$gene}){
			print tarA $gene."\t".$chr{$gene}."\t".$start{$gene}."\t".$end{$gene}."\n";
		    }
		}
		
		if($spegen2{$key1}{$gene} eq $spe2){
		    print tarB $chr{$gene}."\t".$start{$gene}."\t".$end{$gene}."\n";
		}
	    }
	}
	
    }
}
close(tarA);
close(tarB);


my %ono2one=();
my %mutiple=();
my %countChr=();
my $i=0;
my $a=0;
my $k=0;

my $b=0;
my $com=0;

my $spe1uniq=$spe1."_uniq.gff";
my $spe2uniq=$spe2."_uniq.gff";
my $speuniq=$spe1.".".$spe2."_uniq.gff";
open(All,">all_orthologues.gff");
open(Uniq1,">$spe1uniq");
open(Uniq2,">$spe2uniq");
open(Uniq3,">$speuniq");
open(Com,">common.gff");
open(Out,">one2one.txt");
open(Out1,">one2one.loca.txt");
open(Out2,"one2one.gff");

foreach my $key1(sort keys %count){
    #print $key1."\t";
     #foreach my $key2(sort {$a<=>$b}keys $count{$key1}){
     #print $key2."\t";
     #print $count{$key1}{$key2}."\t";
    $b=0;


    foreach my $sp(@species){
	#print "\t$sp";
	#print "$count{$key1}{$sp}\n";
	if($count{$key1}{$sp}>0){
	    $b++;
	    #print "\t$b";
	}
    }
    #print "\n";
    #print $b."\t".scalar(@species)."\n";

    
    my $speAgeneID=$spegen{$key1}{$spe1};
    my $speBgeneID=$spegen{$key1}{$spe2};

    if($chr{$speBgeneID}){
	print All "$chr{$speBgeneID}\tOrthofinder\tgene\t$start{$speBgeneID}\t$end{$speBgeneID}\t.\t$strand{$speBgeneID}\t.\t$speBgeneID\n";
    }

    
    #print $b."\t".scalar(@species)."\n";
    if($b == scalar(@species)-1){
	if($chr{$spegen{$key1}{$spe1}}){
	    print Com "$chr{$speBgeneID}\tOrthofinder\tgene\t$start{$speBgeneID}\t$end{$speBgeneID}\t.\t$strand{$speBgeneID}\t.\t$speBgeneID\n";
	    $com++;
	}
    }
    
    if($b==1){
	if($chr{$spegen{$key1}{$spe2}}){
	    print Uniq2 "$chr{$speBgeneID}\tOrthofinder\tgene\t$start{$speBgeneID}\t$end{$speBgeneID}\t.\t$strand{$speBgeneID}\t.\t$speBgeneID\n";
	}
    }
    if($b==1){
	if($chr{$spegen{$key1}{$spe1}}){
	    print Uniq1 "$chr{$speAgeneID}\tOrthofinder\tgene\t$start{$speAgeneID}\t$end{$speAgeneID}\t.\t$strand{$speAgeneID}\t.\t$speAgeneID\n";
	}
    }

    if($b==2 || $b==1){
	if($chr{$spegen{$key1}{$spe1}} || $chr{$spegen{$key1}{$spe2}}){
	    if($speBgeneID){
		print Uniq3 "$chr{$speBgeneID}\tOrthofinder\tgene\t$start{$speBgeneID}\t$end{$speBgeneID}\t.\t$strand{$speBgeneID}\t.\t$speBgeneID\n";
	    }
	}
    }
    if ($count{$key1}{$spe1} ==1 && $count{$key1}{$spe2} ==1){
	
	if($chr{$spegen{$key1}{$spe1}}){
	    $i++ ;
	    # print $key1." ". $count{$key1}{$spe1}." ".$count{$key1}{$spe2}."\t";
	    # print $i."\n";
	    print Out $chr{$spegen{$key1}{$spe1}}."\t".$spegen{$key1}{$spe1}."\t".$chr{$spegen{$key1}{$spe2}}."\t".$spegen{$key1}{$spe2}."\n";
	    print Out1 $chr{$spegen{$key1}{$spe1}}."\t".$spegen{$key1}{$spe1}."\t".$start{$speAgeneID}."\t".$end{$speAgeneID}."\t".$chr{$spegen{$key1}{$spe2}}."\t".$spegen{$key1}{$spe2}."\t".$start{$speBgeneID}."\t".$end{$speBgeneID}."\n";
	    print Out2 "$chr{$speBgeneID}\tOrthofinder\tgene\t$start{$speBgeneID}\t$end{$speBgeneID}\t.\t$strand{$speBgeneID}\t.\t$speBgeneID\n";
	}
	#print "$chr         \n";
    }else{
	if($count{$key1}{$spe1} >=1 && $count{$key1}{$spe2} >=1){
	    #	print $key1."\n";
	    $k++;
	    #print $k."\n";
	}
    }
	#    }
}

close(Out);
close(Com);
print "$spe1 and $spe2 one to one genefamily number= $i\n";
print "$spe1 and $spe2 Mutiple genefamily number= $k\n";
print "command group number=$com\n";


my %sca2chr=();

open(IN,"one2one.txt")||die "cannot open one2one.txt: $!\n";
while(<IN>){
    #print $_;
    chomp $_;
    @sRows=split/\s+/,$_;
    #print $sRows[2]."\t".$sRows[0]."\n";
    $sca2chr{$sRows[2]}{$sRows[0]}+=1;
}
close(IN);

open(OUT,">final_one2one.txt");
foreach my $k1 (sort keys %sca2chr){
    foreach my $k2 (sort keys %{$sca2chr{$k1}}) {
	if(length($k2)==1){
	    if($k2 != 'Z'){
		print OUT $k1."\t".$k2."\t".$sca2chr{$k1}{$k2}."\n";
	    }else{
		print OUT $k1."\t".$k2."\t".$sca2chr{$k1}{$k2}."\n";
	    }
	}else{
	    print OUT $k1."\t".$k2."\t".$sca2chr{$k1}{$k2}."\n";
	}
    }   
}
close(OUT);

print "Done!\n";


print "outfile1 for heatmap and circos => one2one.txt\tone2one.loca.txt\tfinal_one2one.txt\n";
print "common_link = common.gff\n";
print "All_oligospora = all_orthologues.gff\n";
print "specie1_uniq = $spe1uniq\n";
print "specie2_uniq = $spe2uniq\n";
print "TwoSpecieUniq = $speuniq\n";
