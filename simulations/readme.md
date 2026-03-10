# simulation study

## simulation model

without loss of generality, we can use coalescent units, and aim
to estimate the rate of evolution per coalescent unit.

`L=100` loci

Lévy process for the evolution of each locus effect:
- mutational model (used in the paper):
  * Mutations appear with rate `μ` (Poisson number of mutation),
  * Each mutation has an effect drawn from a centered distribution (mean 0)
    with variance `σ2m`.
    Here we use a centered Laplace distribution
    (aka double exponential distribution):
    `Laplace(0, θ=√(σ2m/2) )` has mean 0 and variance `σ2m`.
- normal (not used in paper, but included in the code): Brownian motion with
  variance rate `σ2=1`

We can choose `μ σ2m = 1`, without loss of generality, that is, `σ2m = 1/μ`.
With `T` being the age of the root, we expect `μT` mutations along each
lineage from the root to a tip (if the phylogeny is time-consistent)
For example, `μT=3` corresponds to an expected 3 mutations along each lineage
between a sampled individual and its ancestor at the root. This setting would
lead to: `μ=3/T` and `σ2m = T/3`, so `θ=√(T/6)` for the Laplace scale parameter.

Then all methods aim to estimate `L μ σ2m = 1`

- variance in the root population: `v0 = μ σ2m λ` with λ ∈ {0, 0.5, 1, 2}
- sample `ni = 20` individuals per population:
  this is representative of many empirical studies
- phylogeny:
  + Polemonium 17-taxon calibrated network `popnet_17taxa.nwk` (used in paper)
  + the 17-population major tree from that network (not used in paper)
- both calibrated tree and network have a total height of 1.
  to vary the amount of ILS, multiply all their branch lengths by `T`,
  with small `T` for more ILS and a large `T` for less ILS.
  T ∈ {1, 10}

100 replicate data sets simulated under each condition.

## code to simulate & analyze data

see `accuracy_GC_simulate.jl`.
needs `popnet_17taxa.nwk` as input for the 17-population species tree
loads `compoundpoisson.jl` for compound Poisson simulation.

## code to visualize the results

see `simulation_figures.qmd`. needs `accuracy_GC_results.csv`,
which is created by `accuracy_GC_simulate.jl`.
