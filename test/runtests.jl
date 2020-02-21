using ClutteredEnvPathOpt

using LightGraphs
using Test

GENERATORS = [cycle_graph, ladder_graph, wheel_graph, x -> grid([x, x])]
RANGE_NODES = 3:32

# @testset "ClutteredEnvPathOpt.jl" begin
#     # Write your own tests here.
# end

include("../src/separator.jl")
@testset "separator tests" begin
    for graph_function in GENERATORS
        for i in RANGE_NODES
            graph = graph_function(i)
            separator = fundamental_cycle_separator(graph, 1)

            # separator is nonempty
            @test length(separator) > 0

            # separator breaks graph into no more than 2 components
            separated_edges = filter(edge -> !(in(src(edge), separator) || in(dst(edge), separator)), collect(edges(graph)))
            separated = SimpleGraphFromIterator(separated_edges)
            components = filter(component -> !(length(component) == 1 && in(component[1], separator)), connected_components(separated))
            @test length(components) < 3
        end
    end
end

@testset "separator postprocessing tests" begin
    for graph_function in GENERATORS
        for i in RANGE_NODES
            graph = graph_function(i)
            separator = fundamental_cycle_separator(graph, 1)
            pp_separator = pp_expell(graph, separator)

            # pp_separator is nonempty
            if length(pp_separator) == 0 println(graph_function, i, separator, pp_separator) end
            @test length(pp_separator) > 0

            # pp_separator breaks graph into no more than 2 components
            separated_edges = filter(edge -> !(in(src(edge), pp_separator) || in(dst(edge), pp_separator)), collect(edges(graph)))
            separated = SimpleGraphFromIterator(separated_edges)
            components = filter(component -> !(length(component) == 1 && in(component[1], pp_separator)), connected_components(separated))
            @test length(components) < 3

            # pp_separator is smaller
            @test length(pp_separator) <= length(separator)
        end
    end
end


@testset "separator helper tests" begin
    @test _find_fundamental_cycle(bfs_parents(cycle_graph(3), 1), 1, Edge(2 => 3)) == Set(1:3)
    @test _find_fundamental_cycle(bfs_parents(cycle_graph(1001), 1), 1, Edge(501 => 502)) == Set(1:1001)
    @test _find_fundamental_cycle(bfs_parents(ladder_graph(8), 1), 1, Edge(15 => 16)) == Set([1, 2, 3, 4, 5, 6, 7, 8, 16, 15])
    
    @test _find_balance(cycle_graph(3), Set(1)) == 0
    @test _find_balance(grid([3, 3]), Set{Int}()) == 0
    @test _find_balance(grid([3, 3]), Set([1, 2, 3])) == 0
    @test _find_balance(grid([3, 3]), Set([4, 5, 6])) == 1
    @test _find_balance(grid([3, 3]), Set([3, 4, 5, 6])) == 2 / 3
    @test _find_balance(ladder_graph(6), Set([2, 8])) == 0.25
end