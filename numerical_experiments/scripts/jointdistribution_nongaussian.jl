#= code to simulate the phenotype of 2 individuals in each population leaf
of a phylogenetic network, conditional on the network's root population P0.
In individuals in P0, the phenotype is assumed to be normally distributed
with mean x0 and variance v0.
=#

using Arrow # to save data frames in smaller files
using CSV
using DataFrames
using Distributions
using PhyloCoalSimulations # requires v1.1.0 or higher
using PhyloNetworks
using PhyloPlots
using Random
using RCall

# for plotting
using AlgebraOfGraphics # v0.11.7. depends on Makie, long installation
# in earlier version: density contours had colors
using CairoMakie # v0.15.6
using PairPlots  # v3.0.3
set_aog_theme!()

numdir = "numerical_experiments"
thisdir = joinpath(numdir, "scripts")
resdir  = joinpath(numdir, "results")

mkdir(resdir)

# 1. setup 3-species network, not time-consistent ---

# set parameters
L = 1 # number of loci
npop = 3
net = readTopology("((b:1,#H1:0):.1,((h:0.1)#H1:0.4::0.5,a:.5):1)P0;")
npop = length(tiplabels(net))
R"pdf"("$resdir/jointdistribution_net_$(npop)pops.pdf", width=3.5, height=2);
R"par"(mar=[0,0,0,0]);
PhyloPlots.plot(net, showedgelength=true, xlim=[.9,4.4], ylim=[.9,3.1],
     style=:majortree, arrowlen=0.1, minorlinetype="longdash", tipoffset=0.1);
R"dev.off"();
transition(x0,len) = Normal(x0, sqrt(len))
# above: evolutionary σ2θ=1 per locus, so σ2L=L×1=1, and from 2N=1
# (lengths in coal units): equilibrium within-species variance = 2N × σ2L =1
# below: choose root variance v0 < equilibrium variance
P0_precision = 10 # 1/v0, called λ in NormalCanon
P0_distribution = NormalCanon(0, P0_precision) # μ = σ²η = 0
nrep = 100_000

# 2. simulate the traits ----------------------------
Random.seed!(342)
X,lab,_ = simulate_polygenictrait(net, nrep, L, P0_distribution, transition;
    nindividuals = Dict("a"=>2, "h"=>2, "b"=>1));
lab == ["b","h_1","h_2","a_1","a_2",] || error("labels in wrong order, check!")
# data frame: 1 row per replicate, 1 column per individual
Xdf = DataFrame(
    a1 = [v[4] for v in X],
    a2 = [v[5] for v in X],
    h1 = [v[2] for v in X],
    h2 = [v[3] for v in X],
    b1 = [v[1] for v in X],
)
# arrow format: 4M vs 9M with CSV.write
# Arrow.write("$resdir/jointdistribution_$(npop)pops.arrow", Xdf; compress=:lz4)
CSV.write(  "$resdir/jointdistribution_$(npop)pops.csv", Xdf)
# Xdf = CSV.read("$resdir/jointdistribution_$(npop)pops.csv", DataFrame)

# 3. create the figure ------------------------------

# to customize, see https://sefffal.github.io/PairPlots.jl/dev/guide
# fg = pairplot(Xdf[!,[:a1,:a2]])

xyl = (; low=-4, high=+4)
xyt = -4:2:4
# fg = Figure(size=(96*5.5, 96*5.5)) # 96px = 1in
fg = Figure(size=(96*7.5, 96*2.5)) # 96px = 1in
pairplot(fg[1,1],
    Xdf[!,[:a2,:a1]] => (
        PairPlots.HexBin(bins=80),
        # PairPlots.Scatter(filtersigma=2, color=:grey60), # to decrease fig file size
        PairPlots.Contour(bandwidth=1.5) ),
    axis=(;
        a1=(; lims=xyl, label="a₁", labelrotation=0),
        a2=(; lims=xyl, label="a₂", labelrotation=0) ),
    bodyaxis=(;
        xticklabelrotation=0, xlabelpadding=-7,
        yticklabelrotation=0, ylabelpadding=-7 )
)
pairplot(fg[1,2], # fg[2,1],
    Xdf[!,[:h2,:h1]] => (
        PairPlots.HexBin(bins=80),
        # PairPlots.Scatter(filtersigma=2, color=:grey60),
        PairPlots.Contour(bandwidth=1.5) ),
    axis=(;
        h1=(; lims=xyl, label="h₁", labelrotation=0),
        h2=(; lims=xyl, label="h₂", labelrotation=0) ),
    bodyaxis=(;
        xticklabelrotation=0, xlabelpadding=-7,
        yticklabelrotation=0, ylabelpadding=-7 )
)
#=
qqnorm(fg[1,2], Xdf[!,:a1],
    markercolor=(:black, 0.3), markersize=1.5,
    axis=(; xticks=xyt, yticks=xyt, ylabelrotation=0,
            xlabel="normal quantiles", ylabel="a₁",)
)
=#
qqnorm(fg[1,3], # fg[2,2],
    Xdf[!,:h1], qqline=:fit, linewidth=0.5, color=:grey,
    markercolor=(:black, 0.3), markersize=1.5,
    axis=(; xticks=xyt, yticks=xyt, ylabelrotation=0,
            xlabel="normal quantiles", ylabel="h₁",)
)
fg

filename = "$resdir/jointdistribution" *
    "_L$(L)_P0prec$(P0_precision)_pairplot.pdf"
save(filename, fg)
