using ClutteredEnvPathOpt

using LightGraphs
using Test


# @testset "ClutteredEnvPathOpt.jl" begin
#     # Write your own tests here.
# end

include("../src/separator.jl")
@testset "separator tests" begin
    small_cycle = cycle_graph(9)
    small_grid = grid([3, 3])
    small_ladder = ladder_graph(5)
    small_wheel = wheel_graph(8)
    big_grid = grid([64, 64])

    pairs = map(
        graph -> (graph, @time fundamental_cycle_separator(graph, 1)),
        [small_cycle, small_grid, small_ladder, small_wheel, big_grid],
    )

    for (graph, separator) in pairs
        # separator is nonempty
        @test length(separator) > 0

        # separator is connected
        is_separator_connected = true
        ordered_vertices = collect(separator)
        for i in 1:length(ordered_vertices)
            for j in 1:length(ordered_vertices)
                if !has_path(graph, ordered_vertices[i], ordered_vertices[j])
                    is_separator_connected = false
                end
            end
        end
        @test is_separator_connected

        # separator breaks graph into no more than 2 components
        separated_edges = filter(edge -> !(in(src(edge), separator) || in(dst(edge), separator)), collect(edges(graph)))
        separated = SimpleGraphFromIterator(separated_edges)
        components = filter(component -> !(length(component) == 1 && in(component[1], separator)), connected_components(separated))
        @test length(components) < 3
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