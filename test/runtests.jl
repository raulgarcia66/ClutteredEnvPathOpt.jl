#ENV["GUROBI_HOME"] = "/Library/gurobi910/mac64"

using ClutteredEnvPathOpt
using LinearAlgebra
using Test
using Pipe
using Plots
using JuMP, Gurobi


# @testset "ClutteredEnvPathOpt.jl" begin
    # Create obstacles
    num_obs = 3;
    seed = 10;
    obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)
    plot();
    ClutteredEnvPathOpt.plot_field(obstacles)
    # points = ClutteredEnvPathOpt.find_points(obstacles)
    # ClutteredEnvPathOpt.plot_points(points)
    display(plot!())
    
    # Set parameters
    N = 30  # number of steps
    f1 = [0.0, 0.1, 0.0]  # initial footstep pose 1
    f2 = [0.0, 0.0, 0.0]  # initial footstep pose 2
    goal = [1.0, 1.0, 0.0]  # goal pose
    Q_g = 10 * Matrix{Float64}(I, 3, 3)  # weight between final footstep and goal pose
    Q_r = [1 0 0; 0 1 0; 0 0 0] # Matrix{Float64}(I, 3, 3)  # weight between footsteps
    q_t = -.05  # weight for trimming unused steps
    # TODO: Assure second footstep is within reachability given first
    
    # Optional named arguments
    d1 = 0.2
    d2 = 0.2
    p1 = [0. 0.07]
    p2 = [0, -0.27]
    delta_x_y_max = 0.1  # max stride norm in space
    delta_θ_max = pi/4  # max difference in θ
    L = 5  # number of pieces in p.w.l approx. of sine/cosine 

    # Compute optimal path
    x1, y1, θ1, t1, stats1 = solve_steps(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="merged")
    x = copy(x1); y = copy(y1); θ = copy(θ1); t = copy(t1);
    term_status, r_solve_time, rel_gap, simplex_iters, r_node_count, num_vertices, merged_size, full_size, num_free_faces, num_free_face_ineq, method_used = stats1
    
    x2, y2, θ2, t2, stats2 = solve_steps(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="full")
    term_status2, r_solve_time2, rel_gap2, simplex_iters2, r_node_count2, num_vertices2, merged_size2, full_size2, num_free_faces2, num_free_face_ineq2, method_used2 = stats2
    
    x3, y3, θ3, t3, stats3 = solve_steps(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="bigM")
    term_status3, r_solve_time3, rel_gap3, simplex_iters3, r_node_count3, num_vertices3, merged_size3, full_size3, num_free_faces3, num_free_face_ineq3, method_used3 = stats3
    
    seed_range = 1:3
    num_obs_range = 1:4
    time_diffs = []
    for seed in seed_range
        for num_obs in num_obs_range
            println("\nOn test Seed = $seed, Num_Obs = $num_obs")
            obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)
            _,_,_,_,stats = solve_steps(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method = "merged")
            _,_,_,_,stats2 = solve_steps(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="full")
            _,_,_,_,stats3 = solve_steps(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="bigM")              
            term_status, solve_time, rel_gap, simplex_iters, node_count, num_vertices, merged_size, full_size, num_free_faces, num_free_face_ineq, method_used = stats
            term_status2, solve_time2, rel_gap2, simplex_iters2, node_count2, num_vertices2, merged_size2, full_size2, num_free_faces2, num_free_face_ineq2, method_used2 = stats2
            term_status3, solve_time3, rel_gap3, simplex_iters3, node_count3, num_vertices3, merged_size3, full_size3, num_free_faces3, num_free_face_ineq3, method_used3 = stats3
            # Merged biclique times vs bigM
            t_diff = solve_time3 - solve_time
            push!(time_diffs, t_diff)
        end
    end

    # Trim excess steps
    num_to_trim = length(filter(tj -> tj > 0.5, t[3:end]))
    if num_to_trim % 2 == 0
        x = vcat(x[1:2], x[num_to_trim + 3 : end]);
        y = vcat(y[1:2], y[num_to_trim + 3 : end]);
        θ = vcat(θ[1:2], θ[num_to_trim + 3 : end]);
    else
        x = vcat(x[1], x[num_to_trim + 3 : end]);
        y = vcat(y[1], y[num_to_trim + 3 : end]);
        θ = vcat(θ[1], θ[num_to_trim + 3 : end]);
    end

    num_to_trim2 = length(filter(tj -> tj > 0.5, t2[3:end]))
    if num_to_trim2 % 2 == 0
        x2 = vcat(x2[1:2], x2[num_to_trim2 + 3 : end]);
        y2 = vcat(y2[1:2], y2[num_to_trim2 + 3 : end]);
        θ2 = vcat(θ2[1:2], θ2[num_to_trim2 + 3 : end]);
    else
        x2 = vcat(x2[1], x2[num_to_trim2 + 3 : end]);
        y2 = vcat(y2[1], y2[num_to_trim2 + 3 : end]);
        θ2 = vcat(θ2[1], θ2[num_to_trim2 + 3 : end]);
    end

    num_to_trim3 = length(filter(tj -> tj > 0.5, t3[3:end]))
    if num_to_trim3 % 2 == 0
        x3 = vcat(x3[1:2], x3[num_to_trim3 + 3 : end]);
        y3 = vcat(y3[1:2], y3[num_to_trim3 + 3 : end]);
        θ3 = vcat(θ3[1:2], θ3[num_to_trim3 + 3 : end]);
    else
        x3 = vcat(x3[1], x3[num_to_trim3 + 3 : end]);
        y3 = vcat(y3[1], y3[num_to_trim3 + 3 : end]);
        θ3 = vcat(θ3[1], θ3[num_to_trim3 + 3 : end]);
    end

    # Plot footstep plan
    plot_steps(obstacles, x, y, θ)
    plot_steps(obstacles, x2, y2, θ2)
    plot_steps(obstacles, x3, y3, θ3)
    png("Path Seed $seed Num Obs $num_obs")


    # Plot intersections of circles
    plot_circles(x, y, θ)
    plot_circles(x2, y2, θ2)
    plot_circles(x3, y3, θ3)

#     @test 1 == 1
# end

# @testset "Optimal biclique comparison" begin
    # obstacles, points, g, faces = ClutteredEnvPathOpt.plot_new(2,"Biclique Comparison Test") # optional seed argument needs to be named
    # skeleton = LabeledGraph(g)

    # tree = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton, faces)
    # (cover_vec, graphviz) = ClutteredEnvPathOpt.tree2digraph(tree)
    # io = open("graphviz.txt", "a")
    # println(io, graphviz)
    # close(io)
    # cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, faces)
    # cover_vec = collect(cover)

#     for i in 1:length(cover_vec)
#         a, b = cover_vec[i]
#         a_vec = collect(a);
#         b_vec = collect(b);

#         x_a = map(n -> points[n].first, a_vec)
#         y_a = map(n -> points[n].second, a_vec)

#         x_b = map(n -> points[n].first, b_vec)
#         y_b = map(n -> points[n].second, b_vec)

#         Plots.scatter(x_a, y_a, color="red", lims=(-0.1, 1.1), series_annotations=(map(n -> Plots.text(string(n), :right, 6, "courier"), a_vec)))
#         Plots.scatter!(x_b, y_b, color="blue", lims=(-0.1, 1.1), series_annotations=(map(n -> Plots.text(string(n), :right, 6, "courier"), b_vec)))

#         # for j in 1:length(obstacles)
#         #     Plots.plot!(obstacles[j], ylim = (0, 1))
#         # end
#         ClutteredEnvPathOpt.plot_field(obstacles)
#         ClutteredEnvPathOpt.plot_lines(obstacles)

#         Plots.savefig(string(i))
#     end

#     feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces))
#     E = @pipe LightGraphs.complement(feg.graph) |> LightGraphs.incidence_matrix(_)
#     lower = ceil(Int, log2(length(faces)))
#     upper = length(cover_vec)

    # plot_optimal(E, obstacles, points, lower, upper)

#     ClutteredEnvPathOpt._is_valid_biclique_cover(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces)), cover)
# end

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

# @testset "biclique cover tests" begin
#     for i in 2:24
#         skeleton = LabeledGraph(ClutteredEnvPathOpt.LightGraphs.grid([i, i]))

#         faces = Set{Vector{Int}}()
#         for j in 1:((i ^ 2) - i)
#             if j % i != 0
#                 face = [
#                     j,
#                     j + 1,      # right
#                     j + i + 1,  # down
#                     j + i,      # left
#                 ]

#                 push!(faces, face)
#             end
#         end

#         cover = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton, faces)

#         @test ClutteredEnvPathOpt._is_valid_biclique_cover(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces)), cover)
#     end

#     for i in 3:24
#         skeleton = LabeledGraph(ClutteredEnvPathOpt.LightGraphs.cycle_graph(i))
#         faces = Set([collect(1:i)])
#         cover = find_biclique_cover(skeleton, faces)

#         @test ClutteredEnvPathOpt._is_valid_biclique_cover(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces)), cover)
#     end
# end
