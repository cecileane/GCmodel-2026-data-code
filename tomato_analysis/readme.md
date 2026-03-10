# tomato data
  * tomato_tree.phy          : time tree from
  [mhibbins/genetreepruningalg](https://github.com/mhibbins/genetreepruningalg/blob/e2627a711e088844f19c680262e05c433ac81a45/analyses/tomato_analysis/make_tomato_test_input_files.py#L13)

  * flower_morphometrics2.csv: floral trait from
  [larabreithaupt/seastaR](https://raw.githubusercontent.com/larabreithaupt/seastaR/refs/heads/main/analyses/flower_morphometrics2.csv)

Tree is re-scaled to coalescent units as in
[larabreithaupt/seastaR](https://raw.githubusercontent.com/larabreithaupt/seastaR/refs/heads/main/analyses/tomato_seastar_analysis.R)

# tomato_analysis_triplets.R

Reproduce the analysis of Hibbins, Breithaupt & Hahn (2023) from
[here](https://raw.githubusercontent.com/larabreithaupt/seastaR/refs/heads/main/analyses/tomato_seastar_analysis.R)

Use phylolm implementation for GC models.

Outputs:
  * results/tomato_high_low.pdf: plot of the triplet analysis (Fig. 8)

# tomato_analysis.R

Analysis of the complete tomato tree, with replicates.

Only one population is kept per species.

Models are fitted using phylolm.

Outputs: 
  * results/tomato_tree_plot.pdf         : plot of the tree and data (Fig. 7)
  * results/phylolm_wholetree_withrep.csv: table with results from the analyses
  * results/tomato_wAIC.pdf              : plot of wAIC weights for all models (Fig. 9)
  * results/tomato_variances.pdf         : plot of estimated evolutionary and within population variances for all models (Fig. 10)
