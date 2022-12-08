use strict;

my $usage='perl AI_index.pl <Metazoan> <nonMetazoan> <GH.hmm> <Orthofinder> <aa>';
# example: perl AI_index.pl $meta $nonmeta $GH Orthogroups.txt APVT.fasta"

# script will use blastdbcmd and NR database to taxify the blast hit, we need to export "ncbi-blast-2.11.0" and "NR location" to environment and put nodesDB.txt to this working folder before running
# export PATH="<ncbi-blast-2.11.0/bin>"
# export PATH="<NR location>"


my $blast1=$ARGV[0];
my $blast2=$ARGV[1];
my $GH_result=$ARGV[2];
my $orthogroup=$ARGV[3];


my %db_ID=();

print "reading GH ...\n";
my $GHhmm=&read_GH($GH_result);
my %GH_Hmm=%$GHhmm;

#print "reading fasta ...\n";
my ($seq,$specie)=&read_fasta($fasta);
my %seq_hash=%$seq;
my %gene2spe=%$specie;


#print "reading Orthofinder ...\n";
my ($OG,$GeneNum)=&read_OG($orthogroup);
my %gene2OG=%$OG;
my %geneNum=%$GeneNum;


#print "reading blast results ...\n";
my($Bit1,$eva1,$hit1,$iden1)=&read_blast($blast1);#metazoan
my %Bit_1=%$Bit1;
my %eva_1=%$eva1;
my %hit_1=%$hit1;
my %iden_1=%$iden1;

my($Bit2,$eva2,$hit2,$iden2)=&read_blast($blast2);#non-metazoan
my %Bit_2=%$Bit2;
my %eva_2=%$eva2;
my %hit_2=%$hit2;
my %iden_2=%$iden2;


#taxified hit
# if taxified output is exist will be skipped
print "taxify blast results ...\n";
my($king1,$order1,$phylum1,$genus1,$specie1)=&hit_taxID($blast1);
my %King_1=%$king1;
my %Order_1=%$order1;
my %Phylum_1=%$phylum1;
my %Genus_1=%$genus1;
my %Specie_1=%$specie1;


my($king2,$order2,$phylum2,$genus2,$specie2)=&hit_taxID($blast2);
my %King_2=%$king2;
my %Order_2=%$order2;
my %Phylum_2=%$phylum2;
my %Genus_2=%$genus2;
my %Specie_2=%$specie2;


my $spe;
if($fasta=~/(\w+)\.fasta/){
    $spe=$1;
}

#estimate AI
open(OUT_H,">$HGT");
#header
print OUT_H "Species\tgeneID\tCAZyme_ID\tOrthogroupID\tDB_type\tAI_index\tHGT_index\tKing\tPhylum\tOrder\tGenus\tSpecie\thitID\tgene_sequence\n";


my $out="$spe"."_AI.txt";
my $HGT="$spe\_HGT\.txt";
my %AI=();
my %HGT=();
my %db_note=();
my $out="$spe"."_AI.txt";

# if gene only have non-metazoan have blast hit, set metazoan evalue to 1
foreach my $spID(sort keys %eva_2){
    if($eva_1{$spID}){
	$AI{$spID}=log($eva_1{$spID}+1e-200)-log($eva_2{$spID}+1e-200);  ## AI index
	$db_note{$spID}="Both";
	$HGT{$spID}=$Bit_2{$spID}-$Bit_1{$spID};   ## bit score difference
    }else{
	$eva_1{$spID}=1;
	$Bit_1{$spID}=0;
	$King_1{$spID}="NA";
	$Order_1{$spID}="NA";
	$Phylum_1{$spID}="NA";
	$Genus_1{$spID}="NA";
	$hit_1{$spID}="NA";
	$AI{$spID}=log(1+1e-200)-log($eva_2{$spID}+1e-200);   ## AI index
	$HGT{$spID}=$Bit_2{$spID}-$Bit_1{$spID};    ## bit score difference
	$db_note{$spID}="OnlyNonMeta";
    }
}

# if gene only have metazoan have blast hit, set non-metazoan evalue to 1
foreach my $spID(sort keys %eva_1){
    if($eva_2{$spID}){
    }else{
	$Bit_2{$spID}=0;
	$AI{$spID}=log($eva_1{$spID}+1e-200)-log(1+1e-200);   ## AI index
	$HGT{$spID}=$Bit_2{$spID}-$Bit_1{$spID};   ## bit score difference
	$db_note{$spID}="OnlyMeta";
	$eva_2{$spID}=1;
	$King_2{$spID}="NA";
        $Order_2{$spID}="NA";
        $Phylum_2{$spID}="NA";
        $Genus_2{$spID}="NA";
	$Specie_2{$spID}="NA";
	$hit_2{$spID}="NA";
    }
}

foreach my $id(sort keys %AI){
    next unless $gene2OG{$id};
    if($egg{$id}){
    }else{
	$egg{$id}="NA";
    }
    if($King_2{$id}){
    }else{
	$King_2{$id}="NA";
    }

    if($Genus_2{$id}){
    }else{
        $Genus_2{$id}="NA";
    }
    
    if($Order_2{$id}){
    }else{
        $Order_2{$id}="NA";
    }

    if($Phylum_2{$id}){
    }else{
        $Phylum_2{$id}="NA";
    }
    
    if($GH_Hmm{$id}){
    }else{
        $GH_Hmm{$id}="not_found";
    }
    my $aa_len=length($seq_hash{$id});
    print OUT_H "$gene2spe{$id}\t$id\t$GH_Hmm{$id}\t$gene2OG{$id}\t$db_note{$id}\t$AI{$id}\t$HGT{$id}\t$King_2{$id}\t$Phylum_2{$id}\t$Order_2{$id}\t$Genus_2{$id}\t$Specie_2{$id}\t$hit_2{$id}\t$seq_hash{$id}\n";
}
close(OUT);
    

print "AI index HGT= $HGT\n";

sub read_blast{
    my $in=$_[0];
    my %Bit=();
    my %evalue=();
    my %hit=();
    my %iden=();
    open(IN,$in)||die "cannot open $in :$!\n";
    while(<IN>){
	next if $_=~/#/;
	chomp $_;
	my @srow=split/\s+/,$_;
	my ($spe,$ID)=split/\|/,$srow[0];
	my $Bit=$srow[11];
	my $eva=$srow[10];
	$Bit{$ID}=$Bit;
	$evalue{$ID}=$eva;
	$hit{$ID}=$srow[1];
	$iden{$ID}=$srow[2];
    }
    close(IN);
    
    return(\%Bit,\%evalue,\%hit,\%iden);
}


#>sp|Q9P7D4|ACON2_SCHPO Homocitrate dehydratase, mitochondrial OS=Schizosaccharomyces pombe (strain 972 / ATCC 24843) OX=284812 GN=SPBP4H10.15 PE=3 SV=3
sub read_fasta{
    my $fasta=$_[0];
    my $seq;
    my $spe;
    my $ID;
    my %seq_hash=();
    my %gene2spe=();
    open(IN,"$fasta")||die "cannot open $fasta:$!\n";
    while(<IN>){
        chomp $_;
        if($_=~/>/){
            if(length($seq)>0){
                $seq_hash{$ID}=$seq;
		$gene2spe{$ID}=$spe;
            }
            my @sRow=split/\s+/,$_;
            $ID=substr($sRow[0],1);
	    ($spe,$ID)=split/\|/,$ID;
	   # print $ID."\n";
            $seq='';
        }else{
            $seq.=$_;

        }
    }
    $seq_hash{$ID}=$seq;
    $gene2spe{$ID}=$spe;
    return(\%seq_hash,\%gene2spe);
    close(IN);
}

sub read_GH{
    my $in=$_[0];
    my %Hmm_GH=();
    open(IN,"$in")||die "cannot open $in :$!\n";
    while(<IN>){
	next if $_=~/\#/;
	my @srow=split/\s+/,$_;
	my $GH=$srow[0];
	substr($GH,-4)='';
	my ($spe,$geneID)=split/\|/,$srow[2];
	$Hmm_GH{$geneID}=$GH;
    }
    return(\%Hmm_GH);
    close(IN);
}


sub read_OG{
    # analyze orthofinder result (orthofinder.txt)
    my $in=$_[0];
    my %OG=();
    my %gene2OG=();
    my %countGene=();
    open(IN,"$in")||die "cannot open $in :$!\n";
    while(<IN>){
	chomp $_;
	my @stmp=split/\s+/,$_;
	my $fam=shift @stmp;
	substr($fam,-1)='';
	foreach my $OG_g(@stmp){
	    my($spe,$gene)=split/\|/,$OG_g;
	    #$OG{$family}{$spe};
	    #$gene2OG{$OG_g}=$fam;
	    $gene2OG{$gene}=$fam;
	    $countGene{$fam}{$spe}+=1;
	}
    }
    return(\%gene2OG,\%countGene);
    close(IN);
}


sub hit_taxID{
    my $in=$_[0];
    my $sub_id=$in.".sub_ids";
    my $command;
    if(-e $sub_id){
    }else{
        $command="cut -f2 $in >$sub_id";
	system($command);
    }
    
    my $tax_map=$sub_id.".map";
    if(-e $tax_map){
    }else{
	## this location may  need to changed as example
	#$command="\/mnt\/nas2\/lihowfun\/bin\/ncbi-blast-2.11.0\/bin\/blastdbcmd -db \/mnt\/nas1\/lihowfun\/db\/nr\/nr -entry_batch $sub_id -outfmt \'\%i \%T\' -out $tax_map";
	$command="blastdbcmd -db nr -entry_batch $sub_id -outfmt \'\%i \%T\' -out $tax_map";
	system($command);
    }
    
    my $taxified=$in.".taxified";
    if(-e $taxified){
	print "$taxified exist\n";
    }else{
	#nodesDB.txt location can be changed as example
        #$command="perl hitID2taxid.pl $in $tax_map \/mnt\/nas2\/lihowfun\/bin\/blobtools\/data\/nodesDB.txt > $taxified";
	$command="perl hitID2taxid.pl $in $tax_map nodesDB.txt > $taxified";
	print "$command\n\n";
	system($command);
    }

    my %gene_tax=();
    my %geneTax_all=();
    my %gene_king=();
    my %gene_order=();
    my %gene_phylum=();
    my %gene_genus=();
    my %gene_specie=();
    open(IN,$taxified)||die "cannot open $taxified :$!\n";
    while(<IN>){
        chomp;
	my($geneID,$taxID,$name_all,$class_all,$king,$order,$phylum,$genus,$specie)=split/\t/,$_;
	$geneTax_all{$geneID}=$name_all;
	$gene_king{$geneID}=$king;
	$gene_order{$geneID}=$order;
	$gene_phylum{$geneID}=$phylum;
	$gene_genus{$geneID}=$genus;
	$gene_specie{$geneID}=$specie;
    }
    return(\%gene_king,\%gene_order,\%gene_phylum,\%gene_genus,\%gene_specie);
}

