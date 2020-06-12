# hare-phylogenomics
Code used for the hare phylogenomics paper where we investigate the occurrence of hybridization throughtout the evolutionary history of hares (genus Lepus).

In each section, I describe the pipeline used in a markdown file (follow links below). Scripts referenced in each section can be found inside the respective folders. 

python scrips use python2.7.

### Sections
1. [Read processing, pseudoreferences and mapping](1.pseudoreferences_and_mapping/1.pseudoreferences_and_mapping.md)
2. [Calling genotypes and creating a consensus fasta](2.call_variants_and_fasta_consensus/2.call_variants_and_fasta_consensus.md)
3. [Species tree and ML whole exome tree](3.species_tree_analysis/3.species_tree_analysis.md)
4. [Divergence time inference (MCMCtree)](4.divergence_time_inference/4.divergence_time_inference.md)
5. [Discordance analyses between gene trees and species tree](5.discordance_analyses/5.discordance_analyses.md)
6. [Network analyses (PhyloNet and Treemix)](6.network_analyses/6.network_analyses.md)
7. [Estimating divergence and diversity per species](7.diversity_divergence_admixture/7.diversity_divergece_admixture.md)
8. [Admixture analyses (Dmin, f-branch, fd, fhom, Dfoil)](8.admixture_analyses/8.admixture_analyses.md)


#### To do's
1. Consensus fasta sequences
2. create_123 python from Divergence time inference with MCMCtree.
3. Organize the MCMC tree inputs/outputs.
5. Remove LGR_SOR from the genomics generals input files.
7. Add remainder of f_hom commands (and regenerate inputs).
