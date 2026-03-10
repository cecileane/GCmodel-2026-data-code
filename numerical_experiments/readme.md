# seastaR_comparison.R

Compare the variance matrices of the `seastaR` model 
with the GC model implemented in `phylolm`
on particular trees.

Outputs: 
  * results/DATE/all_var_A_C_AC.csv                    : triplet variance computation results
  * results/DATE/triplet_caterpillar_comparisons.pdf   : Fig. 3
  * results/DATE/all_var_A_C_AC_AB_CD.csv              : quadruplet variance computation results
  * results/DATE/quadruplet_caterpillar_comparisons.pdf: Fig. 4

# jointdistribution_nongaussian.jl

Explore the joint distribution of the trait at the tips of a triplet,
assuming a single locus with effect evolving like a Brownian Motion.

Outputs: 
  * results/jointdistribution_net_3pops.pdf           : Fig. S1
  * results/jointdistribution_3pops.csv               : data for Fig. S2
  * results/jointdistribution_L1_P0prec10_pairplot.pdf: Fig. S2 
