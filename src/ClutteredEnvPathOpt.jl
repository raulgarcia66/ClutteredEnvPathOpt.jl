module ClutteredEnvPathOpt
using Pipe
import LightGraphs

export LabeledGraph, find_feg_separator_lt, find_separator_lt, find_separator_fcs, find_separator_fcs_best, pp_expell, find_biclique_cover, find_feg_separator_lt_best_root
include("types.jl")
include("separator.jl")
include("biclique.jl")
end # module
