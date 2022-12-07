# Aphelenchoides
## This site provide the codes using in "The Aphelenchoides genomes reveal substantial horizontal gene transfers in the last common ancestor of free-living and major plant parasitic nematodes" figures 1-3


## repeat statistic (fig1 and fig2)
### after running Repeatmasker, we get an <genome fasta name>.out output, then we can use the scripts to count the repeats size  
```
perl repeatmasker.normalized.pl <Repeatmasker output>
perl repeatmasker.stat.pl <normalized output>
```

## Synteny (fig3)
### Synteny relationship between nematodes were inferred using orthofinder, and the scripts below help to generate synteny location
```
perl Orthofinder.one2one_loca.pl Orthogroups.txt SequenceIDs.txt SpeciesIDs.txt <Specie gff A> <Specie gff B>
```
