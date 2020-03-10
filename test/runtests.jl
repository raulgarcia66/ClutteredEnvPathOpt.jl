using ClutteredEnvPathOpt
using Test

@testset "ClutteredEnvPathOpt.jl" begin
    # Write your own tests here.
end

@testset "fundamental cycle separator tests" begin
    GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
    RANGE_NODES = 4:32

    for graph_function in GENERATORS
        for i in RANGE_NODES
            lg = LabeledGraph(graph_function(i))
            (separator, a, b) = find_separator_fcs(lg, 1)

            is_valid = true
            for source in a
                break_early = false
                for destination in b
                    is_path = ClutteredEnvPathOpt.LightGraphs.has_path(lg.graph, source, destination, exclude_vertices=collect(separator))
                    is_valid = is_valid && !is_path
                    if !is_valid
                        break_early = true
                        break
                    end
                end

                if break_early break end
            end
            @test is_valid

            @test all(map(vertex -> in(vertex, union(separator, a, b)), collect(keys(lg.labels))))
        end
    end
end

@testset "separator postprocessing tests" begin
GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
    RANGE_NODES = 4:32

    for graph_function in GENERATORS
        for i in RANGE_NODES
            lg = LabeledGraph(graph_function(i))
            (separator, a, b) = find_separator_fcs(lg, 1)
            (pp_separator, pp_a, pp_b) = pp_expell(lg, separator, a, b)

            is_valid = true
            for source in pp_a
                break_early = false
                for destination in pp_b
                    is_path = ClutteredEnvPathOpt.LightGraphs.has_path(lg.graph, source, destination, exclude_vertices=collect(pp_separator))
                    is_valid = is_valid && !is_path
                    if !is_valid
                        break_early = true
                        break
                    end
                end

                if break_early break end
            end
            @test is_valid

            @test all(map(vertex -> in(vertex, union(pp_separator, pp_a, pp_b)), collect(keys(lg.labels))))

            @test length(pp_separator) <= length(separator)
        end
    end
end

# @testset "separator tests random graphs" begin
#     filenames = [
#         "a280",
#         # "bier127",
#         # "ch130",
#         # "ch150",
#         # "d198",
#     ]
#     for filename in filenames
#         open(("test/delaunay-graphs/$filename.tsp")) do file
#             lines = readlines(file)
#             graph = @pipe match(r"\d+", lines[1]) |> parse(Int, _.match) |> ClutteredEnvPathOpt.LightGraphs.SimpleGraph
    
#             for line in lines[7:end - 1]
#                 edge = map(rm -> parse(Int, rm.match), collect(eachmatch(r"\d+", line)))
#                 ClutteredEnvPathOpt.LightGraphs.add_edge!(graph, edge[2], edge[3]);
#             end
    
#             lg = LabeledGraph(graph)
#             (separator, a, b) = find_separator_fcs(lg, 1)

#             is_valid = true
#             for source in a
#                 break_early = false
#                 for destination in b
#                     is_path = ClutteredEnvPathOpt.LightGraphs.has_path(lg.graph, source, destination, exclude_vertices=collect(separator))
#                     is_valid = is_valid && !is_path
#                     if !is_valid
#                         break_early = true
#                         break
#                     end
#                 end

#                 if break_early break end
#             end
#             @test is_valid
#         end
#     end
# end
