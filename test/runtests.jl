using ClutteredEnvPathOpt
using Test

small_ring = Dict(
    1 => Set([2, 5]),
    2 => Set([1, 3]),
    3 => Set([2, 4]),
    4 => Set([3, 5]),
    5 => Set([4, 1]),
)

small_grid = Dict(
    1 => Set([2, 4]),
    2 => Set([1, 3, 5]),
    3 => Set([2, 6]),
    4 => Set([1, 5, 7]),
    5 => Set([2, 4, 6, 7, 8]),
    6 => Set([3, 5, 9]),
    7 => Set([4, 8]),
    8 => Set([5, 7, 9]),
    9 => Set([6, 8]),
)

cross = Dict(
    1 => Set([2]),
    2 => Set([1, 3, 4, 5]),
    3 => Set([2]),
    4 => Set([2]),
    5 => Set([2]),
)

@testset "ClutteredEnvPathOpt.jl" begin
    # Write your own tests here.
end

@testset "bfs tests" begin
    @test ClutteredEnvPathOpt.bfs(small_ring, 1) == Dict(
        2 => 1,
        3 => 2,
        4 => 5,
        5 => 1,
    )

    @test ClutteredEnvPathOpt.bfs(cross, 1) == Dict(
        2 => 1,
        3 => 2,
        4 => 2,
        5 => 2,
    )

    @test ClutteredEnvPathOpt.bfs(cross, 2) == Dict(
        1 => 2,
        3 => 2,
        4 => 2,
        5 => 2,
    )
end

@testset "edge counting tests" begin
    @test ClutteredEnvPathOpt.edges(small_ring) == Set([
        Set([1, 2]),
        Set([2, 3]),
        Set([3, 4]),
        Set([4, 5]),
        Set([5, 1]),
    ])
    @test ClutteredEnvPathOpt.edges(cross) == Set([
        Set([1, 2]),
        Set([2, 3]),
        Set([2, 4]),
        Set([2, 5]),
    ])
    @test ClutteredEnvPathOpt.tree_edges(Dict(2 => 1, 3 => 2, 4 => 5, 5 => 1)) == Set([
        Set([1, 2]),
        Set([2, 3]),
        Set([4, 5]),
        Set([5, 1]),
    ])
end

@testset "fundamental cycle identifying tests" begin
    @test ClutteredEnvPathOpt.find_fundamental_cycle(ClutteredEnvPathOpt.bfs(small_ring, 1), Set([3, 4])) == keys(small_ring)
end

@testset "find component" begin
    @test ClutteredEnvPathOpt.find_smaller_component(small_ring, Set([1, 2, 3, 4, 5])) == Set()
    @test ClutteredEnvPathOpt.find_smaller_component(small_ring, Set{Int}()) == Set()
    @test ClutteredEnvPathOpt.find_smaller_component(small_grid, Set([3, 4, 5, 6])) == Set([1, 2])
    @test ClutteredEnvPathOpt.find_smaller_component(small_grid, Set([2, 3, 4, 5, 6])) == Set(1)
    @test ClutteredEnvPathOpt.find_smaller_component(small_grid, Set([7, 8, 9])) == Set()

end