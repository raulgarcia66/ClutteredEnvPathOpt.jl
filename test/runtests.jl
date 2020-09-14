using ClutteredEnvPathOpt
using Test
using Pipe

import LightGraphs
import VertexSafeGraphs

@testset "ClutteredEnvPathOpt.jl" begin
    # Write your own tests here.
end

# @testset "Lipton-Tarjan separator on Finite Element Graph test" begin
#     # Grid skeleton
#     for i in 2:32
#         vsg = VertexSafeGraphs.VSafeGraph(ClutteredEnvPathOpt.LightGraphs.grid([i, i]))

#         faces = Set{Set{LightGraphs.Edge}}()
#         for j in 1:((i ^ 2) - i)
#             if j % i != 0
#                 face = Set([
#                     LightGraphs.Edge(j, j + 1),         # left
#                     LightGraphs.Edge(j + 1, j + i + 1), # down
#                     LightGraphs.Edge(j + i + 1, j + i), # right
#                     LightGraphs.Edge(j + i, j),         # up
#                 ])

#                 push!(faces, face)
#             end
#         end

#         (separator, a, b) = find_feg_separator_lt(vsg, faces, 1)

#         @test ClutteredEnvPathOpt._is_valid_separator(vsg, separator, a, b)
#     end

#     # Wheel skeleton
#     for i in 2:32
#         vsg = VertexSafeGraphs.VSafeGraph(ClutteredEnvPathOpt.LightGraphs.wheel_graph(i))

#         faces = Set{Set{LightGraphs.Edge}}()
#         for j in 2:(i - 1)
#             face = Set([
#                 LightGraphs.Edge(1, j),
#                 LightGraphs.Edge(j, j + 1),
#                 LightGraphs.Edge(j + 1, 1),
#             ])

#             push!(faces, face)
#         end

#         (separator, a, b) = find_feg_separator_lt(vsg, faces, 1)

#         @test ClutteredEnvPathOpt._is_valid_separator(vsg, separator, a, b)
#     end
# end

GENERATORS = [
    LightGraphs.cycle_graph,
    LightGraphs.ladder_graph,
    LightGraphs.wheel_graph,
    x -> LightGraphs.grid([x, x])
]

function test_separator_algorithm(algorithm, generators, range)
    for graph_function in generators
        for i in range
            vsg = VertexSafeGraphs.VSafeGraph(graph_function(i))
            (separator, a, b) = algorithm(vsg, 1)

            @test ClutteredEnvPathOpt._is_valid_separator(vsg, separator, a, b)
        end
    end
end

# @testset "Lipton-Tarjan separator tests" begin
#     test_separator_algorithm(find_separator_lt, GENERATORS, 4:10)
# end

# @testset "fundamental cycle separator tests" begin
#     test_separator_algorithm(find_separator_fcs, GENERATORS, 4:32)
# end

# @testset "fundamental cycle best separator tests" begin
#     test_separator_algorithm(find_separator_fcs_best, GENERATORS, 4:10)
# end

# filenames = [
#     "a280",
#     "bier127",
#     "ch130",
#     "ch150",
#     "d198",
#     "d1291",
#     "d1655",
#     "d2103",
#     "d493",
#     "d657",
#     "eil101",
#     "eil51",
#     "eil76",
#     "fl1400",
#     "fl1577",
#     "fl417",
#     "gil262",
#     "kroA100",
#     "kroA150",
#     "kroA200",
#     "kroB100",
#     "kroB150",
#     "kroB200",
#     "kroC100",
#     "kroE100"
# ]

# @testset "fundamental cycle separator tests provided graphs" begin
#     for filename in filenames
#         open(("delaunay-graphs/$filename.tsp.del")) do file
#             lines = readlines(file)
#             graph = @pipe match(r"\d+", lines[1]) |> parse(Int, _.match) |> ClutteredEnvPathOpt.LightGraphs.SimpleGraph
    
#             for line in lines[2:end]
#                 edge = map(rm -> parse(Int, rm.match), collect(eachmatch(r"\d+", line)))
#                 ClutteredEnvPathOpt.LightGraphs.add_edge!(graph, edge[1] + 1, edge[2] + 1); # n + 1 bc graph is 0-indexed
#             end

#             vsg = VertexSafeGraphs.VSafeGraph(graph)
#             (separator, a, b) = find_separator_fcs(vsg, 1)
#             @test ClutteredEnvPathOpt._is_valid_separator(vsg, separator, a, b)
#         end
#     end
# end

# @testset "fundamental cycle separator best tests provided graphs" begin
#     for filename in filenames
#         open(("delaunay-graphs/$filename.tsp.del")) do file
#             lines = readlines(file)
#             graph = @pipe match(r"\d+", lines[1]) |> parse(Int, _.match) |> ClutteredEnvPathOpt.LightGraphs.SimpleGraph
    
#             for line in lines[2:end]
#                 edge = map(rm -> parse(Int, rm.match), collect(eachmatch(r"\d+", line)))
#                 ClutteredEnvPathOpt.LightGraphs.add_edge!(graph, edge[1] + 1, edge[2] + 1); # n + 1 bc graph is 0-indexed
#             end
            
#             vsg = VertexSafeGraphs.VSafeGraph(graph)
#             (separator, a, b) = find_separator_fcs(vsg, 1)
#             @test ClutteredEnvPathOpt._is_valid_separator(vsg, separator, a, b)
#         end
#     end
# end

# @testset "Lipton-Tarjan separator tests provided graphs" begin
#     for filename in filenames
#         open(("delaunay-graphs/$filename.tsp.del")) do file
#             lines = readlines(file)
#             graph = @pipe match(r"\d+", lines[1]) |> parse(Int, _.match) |> ClutteredEnvPathOpt.LightGraphs.SimpleGraph
    
#             for line in lines[2:end]
#                 edge = map(rm -> parse(Int, rm.match), collect(eachmatch(r"\d+", line)))
#                 ClutteredEnvPathOpt.LightGraphs.add_edge!(graph, edge[1] + 1, edge[2] + 1); # n + 1 bc graph is 0-indexed
#             end
            
#             vsg = VertexSafeGraphs.VSafeGraph(graph)
#             (separator, a, b) = find_separator_fcs(vsg, 1)
#             @test ClutteredEnvPathOpt._is_valid_separator(vsg, separator, a, b)
#         end
#     end
# end

# @testset "separator postprocessing tests" begin
# GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
#     RANGE_NODES = 4:32

#     for graph_function in GENERATORS
#         for i in RANGE_NODES
#             vsg = VertexSafeGraphs.VSafeGraph(graph_function(i))
#             (separator, a, b) = find_separator_fcs(vsg, 1)
#             (pp_separator, pp_a, pp_b) = pp_expell(vsg, separator, a, b)

#             @test ClutteredEnvPathOpt._is_valid_separator(vsg, pp_separator, pp_a, pp_b)

#             @test length(pp_separator) <= length(separator)
#         end
#     end
# end

@testset "biclique cover tests" begin
    for i in 3:4
        skeleton = VertexSafeGraphs.VSafeGraph(ClutteredEnvPathOpt.LightGraphs.grid([i, i]))

        faces = Set{Vector{Int}}()
        for j in 1:((i ^ 2) - i)
            if j % i != 0
                face = [
                    j,
                    j + 1,      # right
                    j + i + 1,  # down
                    j + i,      # left
                ]

                push!(faces, face)
            end
        end

        cover = find_biclique_cover(skeleton, faces)
        println(cover)

        @test ClutteredEnvPathOpt._is_valid_biclique_cover(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces)), cover)
    end

    # for i in 3:24
    #     skeleton = LabeledGraph(ClutteredEnvPathOpt.LightGraphs.cycle_graph(i))
    #     faces = Set([collect(1:i)])
    #     cover = find_biclique_cover(skeleton, faces)

    #     @test ClutteredEnvPathOpt._is_valid_biclique_cover(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces)), cover)
    # end
end