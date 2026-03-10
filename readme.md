Code to reproduce simulations and data analysis in
"Trait evolution with incomplete lineage sorting and gene ﬂow: the Gaussian Coalescent model"
(Ané & Bastide, 2026)

Depends on implementations of the GC model in
- the R package [phylolm](https://github.com/lamho86/phylolm) v2.7.0
- and the julia package
  [PhyloTraits](https://github.com/juliaphylo/PhyloTraits.jl) v1.2.0.

Here is an overview of the folders in this repository.
Each folder has a separate `readme.md` for more details.

- `numerical_experiments/scripts/`
  + `seastaR_comparisons.R` to reproduce Figures 3 and 4
  + `jointdistribution_nongaussian.jl` to reproduce Figure S2.

- `simulations/`
  + `accuracy_GC_simulate.jl`: julia code to simulate the data, loading code from
    `compoundpoisson.jl` (for the mutational process) and loading the
    phylogenetic network `popnet_17taxa.nwk` with 17 populations.
  + `simulation_figures.qmd`: R code to visualize the simulation results, and
    reproduce Figures 5-6 and S3-S6.
    This code uses the output csv file created by `accuracy_GC_simulate.jl`.
  + *not* needed to reproduce the simulations: `population_phylogenies.jl`
    to create the input network `popnet_17taxa.nwk` (cleaning that from Teo et al. 2023).
    This script also creates the 45-population tree `poptre_45taxa.nwk`,
    also not needed to reproduce the simulation results.

- `tomato_analysis/` for the analysis of floral traits in wild tomatoes:
  + input tree and trait data in `tomato_tree.phy` and `flower_morphometrics2.csv`
  + `tomato_analysis_triplets.R` to analyze the low and high ILS triplets and
    reproduce Figure 8.
  + `tomato_analysis.R` to analyze the "full" dataset, and reproduce
    Figures 7, 9 and 10.
