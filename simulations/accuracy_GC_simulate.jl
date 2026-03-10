#= simulations to quantify the accuracy to estimate the GC model parameters
1. check accuracy of μ, σ2, v0, and λ when estimated
2. power & type-1 error of LRT to test that λ=1
   and model selection with AIC
=#

dir = "simulations"
cd(dir)

include("compoundpoisson.jl") # loads Distributions & Random
using CSV
using DataFrames
using GLM
using PhyloCoalSimulations # requires v1.1.0
using PhyloNetworks
using PhyloPlots
using PhyloTraits # requires v1.2.0
using RCall
using StatsBase

sigma2_within_gc(m::PhyloNetworkLinearModel, net::HybridNetwork) =
    sigma2_within_gc(lambda_estim(m), sigma2_phylo(m), net)
function sigma2_within_gc(λ, σ2evo, net::HybridNetwork)
    -, eV = PhyloTraits.gaussiancoalescent_covariancematrix(net, λ)
    return σ2evo * eV[:tips]
end

net17 = readnewick("popnet_17taxa.nwk")
# maximum(getnodeheights(net17)) ≈ 1 || error("net17 should have unit height")
tax17 = tiplabels(net17)
ntax = length(tax17) # 17
tre17 = majortree(net17)
# below reticulations: t6 (old), t9 (recent), {t12,t13} (old)
# youngest sister clade: {t3,t4}. outgroup: t17
tax6 = ["t3","t4","t6","t9","t12","t13"]

# 2N=1
topolist = ["net", "tre"] # topology: with / without reticulations
kernlist = ["mutational", "normal"]
Tlist = [10,1]  # multiplier of branch lengths to control ILS
net17_T1  = deepcopy(net17)
net17_T10 = deepcopy(net17)
for e in net17_T10.edge e.length *= 10; end
tre17_T1  = deepcopy(tre17)
tre17_T10 = deepcopy(tre17)
for e in tre17_T10.edge e.length *= 10; end

ni=20 # number of individuals per species
nmutations = 3
L=100
λlist =[0, 0.5, 1, 2]
# σ2m=1
# v0 = λ (2Nσ2m) = λ

#= within-species variation: the implementation matches the formula for Hu:
[[mean(sigma2_within_gc(λ, 1, n)) for λ in λlist] for n in [tre17_T1,tre17_T10]]
[[exp(-T) * λ + (1-exp(-T)) * 1   for λ in λlist] for T in [1,10]]
both output this, up to some precision:
 [0.6321205155486016, 0.8160602577743008, 1.0, 1.3678794844513984]
 [0.9999546000168257, 0.9999773000084126, 1.0, 1.0000453999831744]

but it's higher under a network:
[[mean(sigma2_within_gc(λ, 1, n)) for λ in λlist] for n in [net17_T1,net17_T10]]
 [0.6410849191067046, 0.8322268597952748, 1.0233688004838453, 1.405652681860986]
 [1.0688460659442562, 1.0691905071381087, 1.0695349483319614, 1.0702238307196663]

variation across species with λ=2, T=1:
eV = sigma2_within_gc(2, 1, net17)  # most are = 1.36787944117144
taxa = tiplabels(net17)
[taxa[i] => round(eV[i,1], digits=4) for i in 1:17 if eV[i,1]>1.37]
 "t13" => 1.4826
 "t12" => 1.4826
  "t9" => 1.7227
  "t6" => 1.4257
Its the 4 taxa below some reticulations: all kept in tax6
=#

nrep = 100 # number of replicates

# grab memory to be re-used
df = DataFrame(
    tipnames = Vector{String}(undef, ntax * ni),
    trait    = Vector{Float64}(undef, ntax * ni))
f = @formula(trait ~ 1)
pl_estim = Dict(:lambda => (start=1.0,)) # λ not fixed
pl_fix1  = Dict(:lambda => (start=1.0, fixed=true))
pl_true  = Dict(:lambda => (start=1.0, fixed=true)) # to be modified in loop
ind2pop = r"_\d+$" => "" # removes an ending "_12" for example

nλmodels = 4 # λ estimated, λ fixed at 1, λ fixed at true value, BM+wsv
numrows = length(topolist) * length(Tlist) * length(kernlist) *
          length(λlist) * nrep * nλmodels
res = DataFrame(
    topo   = Vector{String}(undef, numrows),
    netT   = Vector{Float64}(undef, numrows),
    kernel = Vector{String}(undef, numrows),
    lambda = Vector{Float64}(undef, numrows),
    lmodel = Vector{String}(undef, numrows),
    rep    = Vector{Int}(undef, numrows),
    lambda_hat = Vector{Float64}(undef, numrows),
    mu_hat  = Vector{Float64}(undef, numrows),
    s2_hat  = Vector{Float64}(undef, numrows),
    v0_hat  = Vector{Union{Missing,Float64}}(missing, numrows),
    s2_within = Vector{Float64}(undef, numrows),
    loglik  = Vector{Float64}(undef, numrows),
    aic     = Vector{Float64}(undef, numrows),
    lrt_pval = Vector{Union{Missing,Float64}}(missing, numrows),
)
for tax in tax6 # track some populations means (M) and variances (V)
    res[:,Symbol(tax * "_m")] = Vector{Union{Missing,Float64}}(missing, numrows)
    res[:,Symbol(tax * "_v")] = Vector{Union{Missing,Float64}}(missing, numrows)
end

Random.seed!(4184)

# m_gc = nothing # global variable, to look at the last for sanity
irow = 1
for topo in topolist for T in Tlist
  @info "starting topo=$topo, T=$T"
  net = eval(Symbol("$(topo)17_T$T"))
  transition_normal(x,len) = Normal(x, sqrt(len))
  transition = transition_normal # outside scope of next loops
  for kern in kernlist for λ in λlist
    @info "  starting kernel=$kern, λ=$λ"
    pl_true[:lambda] = (start=λ, fixed=true)
    root_prior = Normal(0, sqrt(λ)) # v0 = λ (2Nσ2m) = λ
    μ = nmutations/T
    θ = sqrt(1/(2μ))
    if kern == "mutational"
        # trait variance / time units: σ2m = μ * 2θ² = μ * 1/(μL) = 1
        transition(x,len) = CompoundPoissonLaplace(x, μ*len, θ)
    else  #    same: σ2m = 1
        transition = transition_normal
    end
    x,lab,_ = simulate_polygenictrait(net, nrep, L,
          root_prior, transition; nindividuals=ni);
    for i in eachindex(lab)
        df.tipnames[i] = replace(lab[i], ind2pop)
    end
    for irep in 1:nrep
        df.trait .= x[irep]
        # store populations M & V on row for first λ model
        pMV = combine(groupby(filter(r -> r.tipnames ∈ tax6, df),:tipnames),
                :trait => mean => :m, :trait => var => :v)
        for r in eachrow(pMV)
            tax = r[:tipnames]
            res[irow,Symbol(tax * "_m")] = r[:m]
            res[irow,Symbol(tax * "_v")] = r[:v]
        end
        for λmodel in ("true", "fix1", "estim", "BMw") # order matters for LRT of "λ=1"
            ro = res[irow,:]
            ro[:topo]   = topo
            ro[:netT]   = T
            ro[:kernel] = kern
            ro[:lambda] = λ
            ro[:lmodel] = λmodel
            ro[:rep]    = irep
            estimlam   = (λmodel == "estim")
            model_args = (λmodel == "BMw" ?
                (model = "BM", withinspecies_var = true) :
                (model = "gaussiancoalescent", paramlist = eval(Symbol("pl_$λmodel"))))
            m_gc = phylolm(f, df, net; model_args...)
            ro[:lambda_hat] = lambda_estim(m_gc)
            ro[:mu_hat] = coef(m_gc)[1]
            ro[:s2_hat] = sigma2_phylo(m_gc)
            if λmodel != "BMw"
                ro[:v0_hat] = m_gc.evomodel.v0
            end
            ro[:s2_within] = (λmodel == "BMw" ? sigma2_within(m_gc) :
                mean(sigma2_within_gc(m_gc, net)))
            ro[:loglik] = loglikelihood(m_gc)
            ro[:aic]    = aic(m_gc)
            if estimlam
                # irow-1 = row for "fix1"
                x2 = 2*(ro[:loglik] - res[irow-1,:loglik])
                pval = ccdf(Chisq(1), x2)
                ro[:lrt_pval] = pval
            end
            global irow += 1
        end
    end
  end; end
end; end

CSV.write("accuracy_GC_results.csv", res)
