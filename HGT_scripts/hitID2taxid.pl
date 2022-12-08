use strict;

my $blast=$ARGV[0];
my $taxid=$ARGV[1];
my $nodeDB=$ARGV[2];# names.dmp
my %gene2hit=();
my $hit2tax=&read_taxid($taxid);
my %hit2tax_hash=%$hit2tax;
my($taxName,$taxLV,$tax2up)=&read_nodeDB($nodeDB);
my %taxName_hash=%$taxName;
my %taxLV_hash=%$taxLV;
my %tax2up_hash=%$tax2up;
my %tax2name_all=();
my $class;
my $up;
my $name;
my %tax2_king=();
my $name_all;
my $undef_out=$blast."undef";
my $Eukaryota=$blast."Eukaryota";
my $classided=$blast."classifed_info";


my %tax_class=();
my %tax_name=();
foreach my $tax(sort keys %taxLV_hash){
     $up=$tax2up_hash{$tax};
     $class=$taxLV_hash{$tax};
     #$class_all=$class;
     $name=$taxName_hash{$tax};
     $name_all=$name;
     
     my $oriname=$name;
     my $control="N";
     my $i=0;
     while($control eq "N"){
	 $i++;
	 if($i==1){
	     
	 }else{
	     $up=$tax2up_hash{$up};
	 }

	 $class.="|$taxLV_hash{$up}";
	 #$class_all.=$class;
	 $name=$taxName_hash{$up};
	 $name_all.="|$name";
	 

	 if($class =~/genus/ && $class =~/species/ && $class =~/phylum/ && $class =~/superkingdom/ ){
	     $control="S" ;
         }
	 
	 
	 if($i>20){
	     $control="T";	     
	 }
     }
     if($control eq "S"){
	 $tax_class{$tax}=$class;
	 $tax_name{$tax}=$name_all;
	 $tax2name_all{$tax}=$name_all;
     }elsif($control eq "T"){
	 $tax2_king{$tax}="Too much";
     }elsif($control eq "E"){
	 $tax2_king{$tax}="undef";
     }else{
	 $tax2_king{$tax}="undef";
     }
}


my %genus=();
my %order=();
my %phylum=();
my %superkingdom=();
my %specie=();

my $i=0;
foreach my $tax(sort keys %tax_class){
    my @class_array=split/\|/,$tax_class{$tax};
    my @name_array=split/\|/,$tax_name{$tax};
    $i=0;
    foreach my $class_n(@class_array){
	my $name_now=$name_array[$i];
	$i++;
	next if $name_now eq "";
	if($class_n eq "species"){
	    $specie{$tax}=$name_now;
	}elsif($class_n eq "genus"){
	    $genus{$tax}=$name_now;
	}elsif($class_n eq "family"){
	}elsif($class_n eq "order"){
	    $order{$tax}=$name_now;
	}elsif($class_n eq "phylum"){
	    $phylum{$tax}=$name_now;
	}elsif($class_n eq "superkingdom"){
	    my $b=$i-3;
	    if($name_array[$b] eq "Fungi"){
		$superkingdom{$tax}="Fungi";
	    }else{
		$superkingdom{$tax}=$name_now;
	    }
	}else{
	}
    }
}

open(IN,$blast)||die "cannot open $blast:$!\n";
while(<IN>){
    chomp;
    my($QgeneID,$Sseqid,$pident,$length,$mismatch,$gap,$Qstart,$Qend,$Sstart,$Send,$Eva,$bit)=split/\s+/,$_;
    my $taxID=$hit2tax_hash{$Sseqid};
    my $spe=$taxName_hash{$taxID};
    my $upper=$tax2up_hash{$taxID};
    my $kingdom=$tax2_king{$taxID};
    my $name=$tax2name_all{$taxID};
    my $specie_n=$specie{$taxID};
    my $genus_n=$genus{$taxID};
    my $order_n=$order{$taxID};
    my $phylum_n=$phylum{$taxID};
    my $superkingdom_n=$superkingdom{$taxID};
    my $specie_n=$specie{$taxID};
    my($spe,$geneID)=split/\|/,$QgeneID;
    print "$geneID\t$taxID\t$tax_name{$taxID}\t$tax_class{$taxID}\t$superkingdom_n\t$order_n\t$phylum_n\t$genus_n\t$specie_n\n";
    
}

#1817920	phylum	Candidatus Terrybacteria	1794811
sub read_nodeDB{
    my $in=$_[0];
    my %taxName=();
    my %tax2up=();
    my %taxLV=();
    open(IN,$in)||die "cannot open $in:$!\n";
    while(<IN>){
	chomp $_;
	next if	$_=~/\#/;
        next if $_=~/root/;
	my ($taxid,$level,$name,$up_taxid)=split/\t/,$_;
	$taxName{$taxid}=$name;
	$taxLV{$taxid}=$level;
	$tax2up{$taxid}=$up_taxid;
    }
    return(\%taxName,\%taxLV,\%tax2up);
    close(IN);
}


sub read_taxid{
    my $in=$_[0];
    my %hit2tax=();
    open(IN,$in)||die "cannot open $in:$!\n";
    while(<IN>){
	chomp $_;
	my($hitID,$tax)=split/\s+/,$_;
	my($db,$hit,$jun)=split/\|/,$hitID;
	$hit2tax{$hit}=$tax;
    }
    return(\%hit2tax);
    close(IN);
}

sub read_taxonomy{
    my $in=$_[0];
    my %taxID2Spe=();
    open(IN,"$in")||die "cannot open $in:$!\n";
    while(<IN>){
        next unless $_=~/scientific\ name/;
        my @sRow=split/\|/,$_;
        my $taxID=$sRow[0];
        my $spe=$sRow[1];
        $taxID=~s/^\s+|\s+$//;
        $spe=~s/^\s+|\s+$//;
        $taxID2Spe{$taxID}=$spe;
    }
    return(\%taxID2Spe);
    close(IN);

}
