#= simulations to quantify the accuracy to estimate the GC model parameters,
under a "best-case" scenario:
simulation model almost matching the estimation model

1. check accuracy of μ, σ2, v0, and λ when estimated
2. power & type-1 error of LRT to test that λ=1
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

tretom = readnewick("../tomato_analysis/tomato_tree_trim_coal.phy")
taxtom = tiplabels(tretom)
ntax = length(taxtom) # 12

# 2N=1
topolist = ["tre"] # topology: without reticulations
kernlist = ["mutational"] # mutational kernel only
Ttree = vcv(tretom)[1,1] # =6, total height of the tree
Tlist = [6, 1]  # multiplier of branch lengths to control ILS
tretom_T1  = deepcopy(tretom) # tree re-scaled to unit coal height
for e in tretom_T1.edge e.length /= Ttree; end
tretom_T6 = deepcopy(tretom) # original tree with height 6

ni=3 # number of individuals per species
nmutations = 3
L=100
λlist =[0, 0.5, 1, 2]
# σ2m=1
# v0 = λ (2Nσ2m) = λ

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

Random.seed!(4184)

# m_gc = nothing # global variable, to look at the last for sanity
irow = 1
for topo in topolist for T in Tlist
  @info "starting topo=$topo, T=$T"
  net = eval(Symbol("$(topo)tom_T$T"))
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

CSV.write("accuracy_GC_results_tomato.csv", res)
