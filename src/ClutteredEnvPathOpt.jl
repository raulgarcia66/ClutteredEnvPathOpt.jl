module ClutteredEnvPathOpt
using Pipe
import LightGraphs

export LabeledGraph, find_feg_separator_lt, find_separator_lt, find_separator_fcs, find_separator_fcs_best, pp_expell, find_biclique_cover, find_feg_separator_lt_best
# export solve_deits, plot_steps, plot_circles, plot_points, all other useful functions
include("types.jl")
include("separator.jl")
include("biclique.jl")
include("obstacles.jl")
include("plotting.jl")
include("stepper.jl")
end # module
