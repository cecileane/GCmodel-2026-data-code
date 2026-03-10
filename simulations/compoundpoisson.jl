using Distributions
using Random # for AbstractRNG
import Distributions: rand

"""
    CompoundPoissonLaplace(m, μ, θ)

Compound Poisson distribution with Laplace-distributed effects.
Example:
```jldoctest
julia> using StableRNGs; rng = StableRNG(91);

julia> d = CompoundPoissonLaplace(3, 0.6, 1)

julia> rand(rng, d, 4)
4-element Vector{Float64}:
 3.0
 3.0
 3.0
 4.3636615168889605
```
"""
struct CompoundPoissonLaplace <: Sampleable{Univariate,Continuous}
    "mean"
    m::Float64
    "mean number of events, sampled from Poisson(μ)"
    μ::Float64
    "scale of each event effect, sampled from Laplace(0,θ)"
    θ::Float64
end
function rand(rng::AbstractRNG, s::CompoundPoissonLaplace)
    n = rand(rng, Poisson(s.μ))
    return s.m + sum(rand(rng, Laplace(0,s.θ), n))
end
