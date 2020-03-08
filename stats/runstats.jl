using Plots
using LightGraphs

# Separator statistics
include("../src/separator.jl")
GENERATORS = [x -> LightGraphs.grid([floor(Int, sqrt(x)), floor(Int, sqrt(x))]), x -> ladder_graph(div(x, 2)), cycle_graph, wheel_graph]
RANGE_NODES = 10:500
LABELS = ["Grid" "Ladder" "Cycle" "Wheel"]

data = map(generator -> begin
    graphs = generator.(RANGE_NODES)
    separators = map(graph -> fundamental_cycle_separator(graph, 1), graphs)
    graph_separators = zip(graphs, separators)
    return map(gs -> _find_balance(gs[1], gs[2]), graph_separators)
end, GENERATORS)
plot(RANGE_NODES, hcat(data[1], data[2], data[3], data[4]), title="Balance", label=LABELS)
savefig("b_no_pp.png")
println("Saved b_no_pp")

data = map(generator -> begin
    graphs = generator.(RANGE_NODES)
    pp_separators = map(graph -> pp_expell(graph, fundamental_cycle_separator(graph, 1)), graphs)
    graph_separators = zip(graphs, pp_separators)
    return map(gs -> _find_balance(gs[1], gs[2]), graph_separators)
end, GENERATORS)
plot(RANGE_NODES, hcat(data[1], data[2], data[3], data[4]), title="Balance w/ pp", label=LABELS)
savefig("b_pp.png")
println("Saved b_pp")


data = map(generator -> begin
    graphs = generator.(RANGE_NODES)
    separators = map(graph -> fundamental_cycle_separator(graph, 1), graphs)
    return length.(separators)
end, GENERATORS)
plot(RANGE_NODES, hcat(data[1], data[2], data[3], data[4]), title="Separator Size", label=LABELS)
savefig("s_no_pp.png")
println("Saved s_no_pp")

data = map(generator -> begin
    graphs = generator.(RANGE_NODES)
    pp_separators = map(graph -> pp_expell(graph, fundamental_cycle_separator(graph, 1)), graphs)
    return length.(pp_separators)
end, GENERATORS)
plot(RANGE_NODES, hcat(data[1], data[2], data[3], data[4]), title="Separator Size w/ pp", label=LABELS)
savefig("s_pp.png")
println("Saved s_pp")

function calculate_stats(generator::Function, metric::Function, range_nodes::UnitRange{Int}, title::String, labels::Array{String}, filename::String)
    
end