ENV["GUROBI_HOME"] = "/Library/gurobi910/mac64"

using ClutteredEnvPathOpt
using Test
using Pipe
using JuMP, Gurobi


@testset "ClutteredEnvPathOpt.jl" begin
    # Write your own tests here.
    obstacles, points, g, faces = ClutteredEnvPathOpt.plot_new(2)
    skeleton = LabeledGraph(g)

    tree = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton, faces)
    (cover_vec, graphviz) = ClutteredEnvPathOpt.tree2digraph(tree)
    io = open("graphviz.txt", "a")
    println(io, graphviz)
    close(io)
    # cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, faces)
    # cover_vec = collect(cover)

    for i in 1:length(cover_vec)
        a, b = cover_vec[i]
        a_vec = collect(a);
        b_vec = collect(b);

        x_a = map(n -> points[n].first, a_vec)
        y_a = map(n -> points[n].second, a_vec)

        x_b = map(n -> points[n].first, b_vec)
        y_b = map(n -> points[n].second, b_vec)

        Plots.scatter(x_a, y_a, color="red", lims=(-0.1, 1.1), series_annotations=(map(n -> Plots.text(string(n), :right, 6, "courier"), a_vec)))
        Plots.scatter!(x_b, y_b, color="blue", lims=(-0.1, 1.1), series_annotations=(map(n -> Plots.text(string(n), :right, 6, "courier"), b_vec)))
        
        # for j in 1:length(obstacles)
        #     Plots.plot!(obstacles[j], ylim = (0, 1))
        # end
        ClutteredEnvPathOpt.plot_field(obstacles)
        ClutteredEnvPathOpt.plot_lines(obstacles)

        Plots.savefig(string(i))
    end

    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces))
    E = @pipe LightGraphs.complement(feg.graph) |> LightGraphs.incidence_matrix(_)
    lower = Int(ceil(log2(length(faces))))
    upper = length(cover_vec)

    plot_optimal(E, obstacles, points, lower, upper)

    ClutteredEnvPathOpt._is_valid_biclique_cover(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces)), cover)
end

# @testset "Lipton-Tarjan separator on Finite Element Graph test" begin
#     # Grid skeleton
#     for i in 2:32
#         lg = LabeledGraph(ClutteredEnvPathOpt.LightGraphs.grid([i, i]))

#         faces = Set{Set{Pair{Int, Int}}}()
#         for j in 1:((i ^ 2) - i)
#             if j % i != 0
#                 face = Set([
#                     j => j + 1,         # left
#                     j + 1 => j + i + 1, # down
#                     j + i + 1 => j + i, # right
#                     j + i => j,         # up
#                 ])

#                 push!(faces, face)
#             end
#         end

#         (separator, a, b) = find_feg_separator_lt(lg, faces, 1)

#         @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
#     end

#     # Wheel skeleton
#     for i in 4:64
#         lg = LabeledGraph(ClutteredEnvPathOpt.LightGraphs.wheel_graph(i))

#         faces = Set{Set{Pair{Int, Int}}}()
#         for j in 2:(i - 1)
#             face = Set([
#                 1 => j,
#                 j => j + 1,
#                 j + 1 => 1
#             ])

#             push!(faces, face)
#         end

#         (separator, a, b) = find_feg_separator_lt(lg, faces, 1)

#         @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
#     end
# end

# @testset "Lipton-Tarjan separator tests" begin
#     GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
#     RANGE_NODES = 4:32

#     for graph_function in GENERATORS
#         for i in RANGE_NODES
#             lg = LabeledGraph(graph_function(i))
#             (separator, a, b) = find_separator_lt(lg, 1)

#             @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
#         end
#     end
# end

# @testset "fundamental cycle separator tests" begin
#     GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
#     RANGE_NODES = 4:32

#     for graph_function in GENERATORS
#         for i in RANGE_NODES
#             lg = LabeledGraph(graph_function(i))
#             (separator, a, b) = find_separator_fcs(lg, 1)

#             @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
#         end
#     end
# end

# @testset "fundamental cycle best separator tests" begin
#     GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
#     RANGE_NODES = 4:32

#     for graph_function in GENERATORS
#         for i in RANGE_NODES
#             lg = LabeledGraph(graph_function(i))
#             (separator, a, b) = find_separator_fcs_best(lg, 1)

#             @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
#         end
#     end
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
#             lg = LabeledGraph(graph)
            
#             (separator, a, b) = find_separator_fcs(lg, 1)

#             @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
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
#             lg = LabeledGraph(graph)
            
#             (separator, a, b) = find_separator_fcs_best(lg, 1)

#             @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
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
#             lg = LabeledGraph(graph)
            
#             (separator, a, b) = find_separator_lt(lg, 1)

#             @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
#         end
#     end
# end

# @testset "separator postprocessing tests" begin
# GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
#     RANGE_NODES = 4:32

#     for graph_function in GENERATORS
#         for i in RANGE_NODES
#             lg = LabeledGraph(graph_function(i))
#             (separator, a, b) = find_separator_fcs(lg, 1)
#             (pp_separator, pp_a, pp_b) = pp_expell(lg, separator, a, b)

#             @test ClutteredEnvPathOpt._is_valid_separator(lg, pp_separator, pp_a, pp_b)

#             @test length(pp_separator) <= length(separator)
#         end
#     end
# end

@testset "biclique cover tests" begin
    for i in 2:24
        skeleton = LabeledGraph(ClutteredEnvPathOpt.LightGraphs.grid([i, i]))

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

        cover = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton, faces)

        @test ClutteredEnvPathOpt._is_valid_biclique_cover(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces)), cover)
    end

    for i in 3:24
        skeleton = LabeledGraph(ClutteredEnvPathOpt.LightGraphs.cycle_graph(i))
        faces = Set([collect(1:i)])
        cover = find_biclique_cover(skeleton, faces)

        @test ClutteredEnvPathOpt._is_valid_biclique_cover(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces)), cover)
    end
end