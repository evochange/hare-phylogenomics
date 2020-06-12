# hare-phylogenomics
Pipeline used for analyses of whole-exome data from 15 hare species and 2 outgroup rabbit species. We performed species tree inference, and estimated hybridization with species network inference and summary statistics (D-statistics, fraction of admixture, admixture proportion and f-branch statistics). 

In each section bellow, I describe the pipeline used in a markdown file (follow links below). Scripts referenced in each section can be found inside the respective folders. 

python scrips use python2.7.

**Contact**: mafaldasferreira (at) cibio.up.pt

### Sections
1. [Read processing, pseudoreferences and mapping](1.pseudoreferences_and_mapping/1.pseudoreferences_and_mapping.md)
2. [Calling genotypes and creating a consensus fasta](2.call_variants_and_fasta_consensus/2.call_variants_and_fasta_consensus.md)
3. [Species tree and ML whole exome tree](3.species_tree_analysis/3.species_tree_analysis.md)
4. [Divergence time inference (MCMCtree)](4.divergence_time_inference/4.divergence_time_inference.md)
5. [Discordance analyses between gene trees and species tree](5.discordance_analyses/5.discordance_analyses.md)
6. [Network analyses (PhyloNet and Treemix)](6.network_analyses/6.network_analyses.md)
7. [Estimating divergence and diversity per species](7.diversity_divergence_admixture/7.diversity_divergece_admixture.md)
8. [Admixture analyses (_D<sub>min</sub>_, f-branch, _f<sub>d</sub>_, _f<sub>hom</sub>_)](8.admixture_analyses/8.admixture_analyses.md)

### Links to the software used in this pipeline 

*Note*: This is not an exhaustive list and some may be missing. All software should be also detailed in each section.

#### General data processing
- [expHTS](https://github.com/msettles/expHTS)
- [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)
- [samtools and bcftools](http://www.htslib.org)
- [Picard Tools](https://broadinstitute.github.io/picard/)
- [Genome Analysis Toolkit](https://gatk.broadinstitute.org/hc/en-us)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)

#### Pseudo-reference generation
- [pseudo-it](https://github.com/bricesarver/pseudo-it)

#### Alignment processing
- [AMAS](https://github.com/marekborowiec/AMAS)
- [TriFusion](https://github.com/ODiogoSilva/TriFusion)
- R package [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
- [snp-sites](https://github.com/sanger-pathogens/snp-sites)
- msa-split, from [phast](http://compgen.cshl.edu/phast/)
- [Genomics General](https://github.com/simonhmartin/genomics_general)

#### Phylogenomic or Population Genomic analyses
- [ASTRAL-III](https://github.com/smirarab/ASTRAL/)
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)
- [PAUP*](https://paup.phylosolutions.com) (implements SVDquartets)
- [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) (implements MCMCtree)
- [SplitsTree4](http://www.splitstree.org)
- [PhyloNet](https://bioinfocs.rice.edu/phylonet)
- [TreeMix](https://bitbucket.org/nygcresearch/treemix/wiki/Home)
- [Genomics General](https://github.com/simonhmartin/genomics_general)
- R package [ape](https://cran.r-project.org/web/packages/ape/index.html)
- R package [phangorn](https://github.com/KlausVigo/phangorn)
- R package [treeman](https://github.com/DomBennett/treeman)


