# Aphelenchoides
## This site provide the codes for the genome project of Aphelenchoides nematodes (The Aphelenchoides genomes reveal substantial horizontal gene transfers in the last common ancestor of free-living and major plant parasitic nematodes).

&nbsp;
&nbsp;

## Statistic of repeat elements (fig1 and fig2)
### after running Repeatmoduler and TransposonPSI we combined the fasta file together to cluster the consensus repeat elements using USEARCH. And the consensus of repeat fasta files were used to compute the repeat locations using Repeatmasker. The output of Repeatmasker <genome fasta name>.out can use the bedtools merge tool to merge the overlap regions and then used the scrips below to count the repeats size.  

&nbsp;
### RepeatModuler, TransposonPSI, USEARCH and Repeatmasker commonds are following:  
```
# RepeatModuler 
RepeatModeler -engine ncbi -pa 5 -database <genome fasta Library>  # use BuildDatabase to creat genome fasta Library first

# TransposonPSI
transposonPSI.pl <genome fasta> nuc
  
# USEARCH
usearch8.1.1861_i86linux64 -sortbylength merged.fa -fastaout seqs_sorted.fasta # merged.fa is the file cat from RepeatModuler and TransposonPSI output  
usearch8.1.1861_i86linux64 -strand both -cluster_smallmem seqs_sorted.fasta -id 0.8 -centroids nr.fasta -uc clusters.uc  
  
# Repeatmasker    
mkdir denovoRepeat  
RepeatMasker -pa <CPU number> -gff  -dir denovoRepeat  -lib repeatLib.fa  -xsmall <genome fasta> # repeatLib.fa is the repeat library created by the USEARCH  
```  
### Statistic  
```
perl repeatmasker.normalized.pl <Repeatmasker output>
perl repeatmasker.stat.pl <normalized output>
```
&nbsp;
&nbsp;
## Batch running and Count the copy numbers of CAZyme (fig2)
### This script help to run the CAZyme script of mutiple species in one script. We put all the protein fasta file in same folder in <specie_name>.fasta format. After running this script, the copies number of CAZyme across each nematodes species with output name "HmmNum.txt" will be generated.

```
perl runCAZYme.pl . 
```
&nbsp;
&nbsp;
## Synteny (fig3)
### Synteny relationship between nematodes were inferred using orthofinder, and the scripts below help to generate synteny location. The input files from orthofinder included Orthogroups.txt SequenceIDs.txt SpeciesIDs.txt. 
#### Code might need to be modified according to the format of gff

```
perl Orthofinder.one2one_loca.pl Orthogroups.txt SequenceIDs.txt SpeciesIDs.txt <Specie gff A> <Specie gff B>
``` 
&nbsp;
&nbsp;
## Horizontal gene transfer (HGT) (fig4 and fig5)
### To estimate the HGT possibility of each genes across nematodes, we applied Alienness Index (AI) (Rancurel et al., 2017; http://alienness.sophia.inra.fr/cgi/index.cgi) method. To speed up the running process, we sepatated the Metazoan and non-Metazoan into two fasta files from NCBI NR database using NCBI tool (ncbi-blast-2.11.0).
&nbsp;
ncbi-blast-2.11.0 is need to install and export (ncbi-blast-2.11.0/bin) to environment  
#### The example of Metazoan database were generated following the command below:  
  
```  
get_species_taxids.sh -n metazoans # search the taxID of Metazoan
get_species_taxids.sh -t 33208 > metazoans.33208.txids # get all the Metazoan taxID into file, you may need to exclude the taxID which are closely related with the species we used to prevent self-alignment
blastdbcmd -db <NR database location> -taxidlist metazoans.33208.txids -target_only -out metazoan.fa  # get the metazoan output according to the taxID
diamond makedb -d metazoan --in metazoan.fa # diamond database format  
```
&nbsp;  
### AI score of genes across nematodes
#### To combine the CAZyme results and AI index together, we worked in the same folder of CAZyme (runCAZYme.pl). Before running AI_index, we need to set up the environment and prepare the input data. We used the ncbi tool (blastdbcmd) and assigned the hit phyla according to the nodeDB.txt (generate from: https://github.com/blaxterlab/blobology)
 
#### The following script and file need to put in this folder or export to environment   
* blast_meta.diamond.pl
* hitID2taxid.pl
* AI_index.pl  
* AI_index.batch.pl  
&nbsp;   
* Orthofinder.txt (from orthofinder)
* metazoan and non-dmnd (from diamond makedb)
* ncbi-blast-2.11.0/bin 
* nodesDB.txt  (https://drive.google.com/file/d/11HHqQslWfJ7Id8fqV8hLJ2DDFFRkO_pR/view?usp=share_link; we generated in the beginning of 2022, it should be updated when NR have better version)
  
 ```
export PATH="HGT_scrpts" # export the environment HGT scripts  
  perl AI_index.batch.pl . <Metazoan.dmnd> <non-Metazoan.dmnd>
  
 ```
  
 &nbsp;  
 output format of AI information will be named <Species name>_HGT according to the format of protein fasta (<Species name>.fasta), and it can be merged into one file
 ```
cat *_HGT > combind_HGT.txt  
 ```
  
  
  
  
  
  

