# Aphelenchoides
## This site provide the codes for the genome project of Aphelenchoides nematodes (The Aphelenchoides genomes reveal substantial horizontal gene transfers in the last common ancestor of free-living and major plant parasitic nematodes)


## Statistic of repeat elements (fig1 and fig2)
### after running Repeatmoduler and TransposonPSI we combined the fasta file together to cluster the consensus repeat elements using USEARCH. And the consensus of repeat fasta files were used to compute the repeat locations using Repeatmasker. The output of Repeatmasker <genome fasta name>.out can use the bedtools merge tool to merge the overlap regions and then used the scrips below to count the repeats size  
```
perl repeatmasker.normalized.pl <Repeatmasker output>
perl repeatmasker.stat.pl <normalized output>
```

## Synteny (fig3)
### Synteny relationship between nematodes were inferred using orthofinder, and the scripts below help to generate synteny location. The input files from orthofinder included Orthogroups.txt SequenceIDs.txt SpeciesIDs.txt. Code might need to modify according to the format of gff

```
perl Orthofinder.one2one_loca.pl Orthogroups.txt SequenceIDs.txt SpeciesIDs.txt <Specie gff A> <Specie gff B>
``` 
  
## Horizontal gene transfer (fig4 and fig5)
### To estimate the genome  
