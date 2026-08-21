# seastaR_comparison.R

Compare the variance matrices of the `seastaR` model 
with the GC model implemented in `phylolm`
on particular trees.

Outputs: 
  * results/DATE/all_var_A_C_AC.csv                    : triplet variance computation results
  * results/DATE/triplet_caterpillar_comparisons.pdf   : Fig. 3
  * results/DATE/all_var_A_C_AC_AB_CD.csv              : quadruplet variance computation results
  * results/DATE/quadruplet_caterpillar_comparisons.pdf: Fig. 4

Further outputs when associated with script
`seastaR_comparisons_sample_gene_trees.jl` (see below):
  * results/DATE/trees/*                                       : Species trees and simulated gene trees 
  * results/DATE/triplet_caterpillar_comparisons_samples.pdf   : Fig. S3
  * results/DATE/quadruplet_caterpillar_comparisons_samples.pdf: Fig. S4
  * results/DATE/all_caterpillar_comparisons_time.pdf          : Fig. S5

# seastaR_comparisons_sample_gene_trees.jl

Julia script to simulate gene trees under the network multispecies coalescent,
along the species trees in `numerical_experiments/results/DATE/trees`,
produced by `seastaR_comparisons.R` (in this folder).
The simulated gene trees are written to
`results/DATE/trees/{triplet,quadruplot}_gene_trees_*`.

# jointdistribution_nongaussian.jl

Explore the joint distribution of the trait at the tips of a triplet,
assuming a single locus with effect evolving like a Brownian Motion.

Outputs: 
  * results/jointdistribution_net_3pops.pdf           : Fig. S1
  * results/jointdistribution_3pops.csv               : data for Fig. S2
  * results/jointdistribution_L1_P0prec10_pairplot.pdf: Fig. S2 
