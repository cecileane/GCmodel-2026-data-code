#= Simulation of gene trees under the network multispecies coalescent.

Input: species trees produced by R script `numerical_experiments/scripts/seastaR_comparisons.R`
* results/DATE/trees/triplet_tree_*   : trees from example 1 (nested caterpillar trees - three tips)
* results/DATE/trees/quadruplet_tree_*: trees from example 2 (mixed setting with ladder - four tips)

Output: simulated gene trees
* results/DATE/trees/triplet_gene_trees_*
* results/DATE/trees/quadruplet_gene_trees_*

=#

dir = "numerical_experiments/results/2026-07-28/trees"
cd(dir)

using PhyloCoalSimulations # requires v1.1.0
using PhyloNetworks
using Random

ngenetrees = 1000

#= Example 1: nested caterpillar tree - three tips
This loop can only be run AFTER lines 120-135 of R script 
`numerical_experiments/scripts/seastaR_comparisons.R`
have been executed.
=#

Random.seed!(1661)
for Ttot in [3, 30]
    for n in [3, 10, 20]
        for t in 0:0.1:1
            tree = readnewick("triplet_tree_n_$(n)_t_$(round(t*Ttot, digits=2))_Ttot_$(Ttot).tree")
            gene_trees = simulatecoalescent(tree, ngenetrees, 1)
            writemultinewick(gene_trees, "triplet_gene_trees_n_$(n)_t_$(round(t*Ttot, digits=2))_Ttot_$(Ttot).tree")
        end
    end
end

#= Example 2: mixed setting with ladder - four tips
This loop can only be run AFTER lines 520-536 of R script
`numerical_experiments/scripts/seastaR_comparisons.R`
have been executed.
=#
Random.seed!(1715)
for Ttot in [3, 30]
    for n in [4, 10, 20]
        for t in 0:0.1:0.9
            tree = readnewick("quadruplet_tree_n_$(n)_t_$(round(t*Ttot, digits=2))_Ttot_$(Ttot).tree")
            gene_trees = simulatecoalescent(tree, ngenetrees, 1)
            writemultinewick(gene_trees, "quadruplet_gene_trees_n_$(n)_t_$(round(t*Ttot, digits=2))_Ttot_$(Ttot).tree")
        end
    end
end
