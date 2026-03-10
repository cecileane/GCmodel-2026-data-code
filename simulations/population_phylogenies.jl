#=
network and tree used for simulations:
from the flowering plant genus Polemonium,
calibrated to fit average genetic distances from gene trees,
with a total height (from root to tips) of 1,
by Teo et al. 2023: https://doi.org/10.18061/bssb.v1i8.8977

`45taxa_caltree.tre` and `17taxa_calnet.tre`
available on dryad: https://doi.org/10.5061/dryad.9ghx3ffkc

This script renames species (and morphs) to:
t1 ... t17 for those in both the tree and the network, and
t18 ... t45 for those in the tree only.
Edge lengths are rounded to 6 decimal places.
In the tree, remove node labels, which were from ASTRAL's bootstrap support

The resulting phylogenies are saved, in (extended) newick format in
`popnet_17taxa.nwk` and `poptre_45taxa.nwk`.
=#
using PhyloNetworks
using PhyloTraits

dir = "simulations"
cd(dir)
tre45 = readnewick("45taxa_caltree.tre")
net17 = readnewick("17taxa_calnet.tre")

lab45 = tiplabels(tre45)
lab17 = tiplabels(net17)
for i in [-10,-12,-13,-16,-20,-22] rotate!(net17, i); end
for i in [47,86,75,68,80,82,59] rotate!(tre45, i); end

lab17 ⊆ lab45 || error("the network should have a subsample of taxa")
# d = Dict(lab => "t$(18-i)" for (i,lab) in enumerate(lab17))
d = Dict(
  "chartaceum_2"             => "t1",
  "eddyense"                 => "t2",
  "brandegeei"               => "t3",
  "confertum"                => "t4",
  "foliosissimum_2"          => "t5",
  "apachianum"               => "t6",
  "pauciflorum"              => "t7",
  "aff._viscosum_sp._nov._2" => "t8",
  "elusum_2"                 => "t9",
  "carneum_2"                => "t10",
  "elegans"                  => "t11",
  "pulcherrimum_shastense"   => "t12",
  "delicatum_2"              => "t13",
  "pectinatum"               => "t14",
  "reptans_2"                => "t15",
  "occidentale_occidentale"  => "t16",
  "micranthum"               => "t17",
  # below: taxa in tree only
  "eximium_2" => "t18",
  "chartaceum" => "t19",
  "eximium" => "t20",
  "viscosum" => "t21",
  "foliosissimum_3" => "t22",
  "foliosissimum" => "t23",
  "albiflorum" => "t24",
  "filicinum" => "t25",
  "flavum" => "t26",
  "aff._viscosum_sp._nov." => "t27",
  "elusum" => "t28",
  "nevadense" => "t29",
  "carneum" => "t30",
)
for n in net17.node
    n.leaf || continue
    n.name = d[n.name]
end
writenewick(net17, "popnet_17taxa.nwk"; round=true, digits=6)

i=45
for n in tre45.node
    if n.leaf
        if !haskey(d, n.name)
            d[n.name] = "t$i"
            i -= 1
        end
        n.name = d[n.name]
    elseif !n.hybrid
        n.name = ""
    end
end
writenewick(tre45, "poptre_45taxa.nwk"; round=true, digits=6)
