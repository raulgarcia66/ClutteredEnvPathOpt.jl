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
    x1, y1, θ1, t1, stats1 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="merged")
    x = copy(x1); y = copy(y1); θ = copy(θ1); t = copy(t1);
    term_status, r_solve_time, rel_gap, simplex_iters, r_node_count, num_vertices, merged_size, full_size, num_free_faces, num_free_face_ineq, method_used = stats1
    
    x2, y2, θ2, t2, stats2 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="full")
    term_status2, r_solve_time2, rel_gap2, simplex_iters2, r_node_count2, num_vertices2, merged_size2, full_size2, num_free_faces2, num_free_face_ineq2, method_used2 = stats2
    
    x3, y3, θ3, t3, stats3 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="bigM")
    term_status3, r_solve_time3, rel_gap3, simplex_iters3, r_node_count3, num_vertices3, merged_size3, full_size3, num_free_faces3, num_free_face_ineq3, method_used3 = stats3
    
    seed_range = 1:3
    num_obs_range = 1:4
    time_diffs = []
    for seed in seed_range
        for num_obs in num_obs_range
            println("\nOn test Seed = $seed, Num_Obs = $num_obs")
            obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)
            _,_,_,_,stats = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method = "merged")
            _,_,_,_,stats2 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="full")
            _,_,_,_,stats3 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="bigM")              
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


####################################### Unofficial Tests ########################################

# Problem Data
num_obs = 2;
seed = 3;
obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs, "Obstacles Seed $seed Num Obs $num_obs", seed = seed)
skeleton = LabeledGraph(g)
all_faces = union(obstacle_faces, free_faces)

# Plot obstacles
plot(title="Obstacles");
ClutteredEnvPathOpt.plot_field(obstacles)
ClutteredEnvPathOpt.plot_lines(obstacles)
ClutteredEnvPathOpt.plot_borders()
ClutteredEnvPathOpt.plot_intersections(obstacles)
png("Obstacles Seed $seed Num Obs $num_obs")
display(plot!())

# Plot faces
ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, plot_name = "Obstacle Faces", col = "dodgerblue")
png("Obstacle Faces Seed $seed Num Obs $num_obs")

ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "Free Faces")
png("Free Faces Seed $seed Num Obs $num_obs")

ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, plot_name = "All Faces", col = "red")
ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "All Faces", col = "green", new_plot = false)
png("All Faces Seed $seed Num Obs $num_obs")

# Plot faces individually
ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, col = "dodgerblue", individually = true)
ClutteredEnvPathOpt.plot_faces(free_faces, points, individually = true)

# Plot edges
ClutteredEnvPathOpt.plot_edges(skeleton, points, obstacles, plot_name = "Skeleton")
png("Skeleton Seed $seed Num Obs $num_obs")

# feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces))
# ClutteredEnvPathOpt.plot_edges(feg, points, obstacles, plot_name = "Finite Element Graph")
# png("FEG Seed $seed Num Obs $num_obs")

feg_S = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
ClutteredEnvPathOpt.plot_edges(feg_S, points, obstacles, plot_name = "Finite Element Graph of S",col="black")
png("FEG of S Seed $seed Num Obs $num_obs")

# feg_obs = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(obstacle_faces))
# ClutteredEnvPathOpt.plot_edges(feg_obs, points, obstacles, plot_name = "Finite Element Graph Obstacle Faces")
# png("FEG Obstacles Seed $seed Num Obs $num_obs")


# Duplicate tests
# obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)
# points, mapped, inside_quant = ClutteredEnvPathOpt.find_intersections(obstacles)
# neighbors = ClutteredEnvPathOpt._find_neighbors(points, mapped)

# dup, dup_ind = face_duplicates(obstacle_faces)
# dup, dup_ind = face_duplicates(free_faces)
# dup, dup_ind = face_duplicates(all_faces)
# plot_intersections_indiv(points)
# dup, dup_ind = point_duplicates(points)
# dup, dup_ind = points_in_hs_duplicates(mapped)
# dup, dup_ind = hs_duplicates(mapped)
# points_in_mapped_ordered(obstacles, mapped)
# plot_neighbors(obstacles, neighbors, points)
# list_neighbors(neighbors, points)
# suspect = suspect_neighbors(neighbors)



# cover3 = ClutteredEnvPathOpt._find_biclique_cover_visual(skeleton, free_faces, points)
# tree = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton, free_faces)

# Free Faces Cover
cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces)
(valid_cover, size_diff, missing_edges, missing_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces)), cover)

feg_free = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg_free)
(valid_cover, size_diff, missing_edges, missing_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces)), merged_cover)

ClutteredEnvPathOpt.plot_edges(missing_edges_lg, points, obstacles, plot_name = "Missing Edges")
png("Missing Edges Seed $seed Num Obs $num_obs Free")

file_name = "Biclique Cover Debug Output Seed $seed Num Obs $num_obs Free Unedited.txt"
f = open(file_name, "w")
write(f, "Plot 1\n")
flush(f)
close(f)
cover2 = ClutteredEnvPathOpt.find_biclique_cover_debug(skeleton, free_faces, points, obstacles, file_name)

# To produce individual FEGs/Skeletons in the recursion process
(C,A,B), skeleton_ac, faces_ac, skeleton_bc, faces_bc = ClutteredEnvPathOpt.find_biclique_cover_one_iter(skeleton, free_faces)
feg_S_ac = ClutteredEnvPathOpt._find_finite_element_graph(skeleton_ac, ClutteredEnvPathOpt._find_face_pairs(faces_ac))
ClutteredEnvPathOpt.plot_edges(feg_S_ac, points, obstacles, plot_name = "Finite Element Graph of S_ac",vertices=feg_S_ac.labels,col="black")
# plot!(title="",axis=([], false))
png("Ex FEG of S_ac Seed $seed Num Obs $num_obs")
feg_S_bc = ClutteredEnvPathOpt._find_finite_element_graph(skeleton_bc, ClutteredEnvPathOpt._find_face_pairs(faces_bc))
ClutteredEnvPathOpt.plot_edges(feg_S_bc, points, obstacles, plot_name = "Finite Element Graph of S_bc", vertices=feg_S_bc.labels, col="black")
# plot!(title="",axis=([], false))
png("Ex FEG of S_bc Seed $seed Num Obs $num_obs")

# DON'T NEED TO ALL FACES COVER; JUST GIVING THE FREE FACES IS VALID
# All Faces Cover
# cover3 = ClutteredEnvPathOpt.find_biclique_cover(skeleton, all_faces)
# (valid_cover, size_diff, missing_edges, missing_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces)), cover)
# ClutteredEnvPathOpt.plot_edges(missing_edges_lg, points, obstacles, plot_name = "Missing Edges")
# png("Missing Edges Seed $seed Num Obs $num_obs All")

# file_name = "Biclique Cover Debug Output Seed $seed Num Obs $num_obs All Unedited.txt"
# f = open(file_name, "w")
# write(f, "Plot 1\n")
# flush(f)
# close(f)
# cover4 = ClutteredEnvPathOpt.find_biclique_cover_debug(skeleton, all_faces, points, obstacles, file_name)


# Biclique Cover Validity Tests-----------------------------------------------------------------------
seed_range = 1:100
num_obs_range = 1:4
fails, tuples = biclique_cover_validity_tests(seed_range, num_obs_range, faces = "all")
fails, tuples = biclique_cover_validity_tests(seed_range, num_obs_range, faces = "free")

percent_decreases = biclique_cover_merger_stats(seed_range, num_obs_range)


# Biclique validity tester function
function biclique_cover_validity_tests(seed_range, num_obs_range; faces::String, with_plots=false)
    fails = 0
    tuples = []
    for seed = seed_range
        for num_obs = num_obs_range
            println("On test Seed = $seed, Num_Obs = $num_obs")
            obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs,"Obstacles Seed $seed Num Obs $num_obs",seed=seed); # optional seed argument, default is 11
            skeleton = LabeledGraph(g)
            all_faces = union(obstacle_faces, free_faces)
            if faces == "all"
                cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, all_faces);
                feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces));
                merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
                # (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
                (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover);
                println("Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")
            else
                cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces);
                feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces));
                merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
                # (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
                (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover);
                println("Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")
            end
            if size_diff != 0
                fails += 1
                push!(tuples, (seed, num_obs))
                println("Test Failed")

                if with_plots
                    plot(title="Obstacles");
                    ClutteredEnvPathOpt.plot_field(obstacles)
                    ClutteredEnvPathOpt.plot_lines(obstacles)
                    ClutteredEnvPathOpt.plot_borders()
                    ClutteredEnvPathOpt.plot_intersections(obstacles);
                    display(plot!())
                    png("Obstacles Seed $seed Num Obs $num_obs")

                    ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "Free Faces")
                    png("Free Faces Seed $seed Num Obs $num_obs")

                    ClutteredEnvPathOpt.plot_edges(skeleton, points, obstacles, plot_name = "Skeleton")
                    png("Skeleton Seed $seed Num Obs $num_obs")

                    if size_diff > 0
                        ClutteredEnvPathOpt.plot_edges(miss_edges_lg, points, obstacles, plot_name = "Missing Edges")
                        if faces == "all"
                            png("Missing Edges All Faces Used Seed $seed Num Obs $num_obs")
                        else
                            png("Missing Edges Free Faces Used Seed $seed Num Obs $num_obs")
                        end
                    end
                end
            end
        end
    end
    return (fails, tuples)
end

# Biclique cover merger function for stats
function biclique_cover_merger_stats(seed_range, num_obs_range)
    file_name = "Biclique Cover Merging Stats.txt" #Seed Range = $seed_range, Num Obs Range = $num_obs_range.txt"
    f = open(file_name, "w")
    write(f, "Seed Range = $seed_range, Num Obs Range = $num_obs_range\n")
    flush(f)
    close(f)
    
    f = open(file_name, "a")
    fails = 0
    tuples = []
    percent_vec = []
    for seed = seed_range
        for num_obs = num_obs_range
            # f = open(file_name, "a")
            write(f, "\nTest Seed = $seed, Num_Obs = $num_obs\n")
            # flush(f)
            # close(f)

            # println("On test Seed = $seed, Num_Obs = $num_obs")
            obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs,"Obstacles Seed $seed Num Obs $num_obs",seed=seed); # optional seed argument, default is 11
            skeleton = LabeledGraph(g)
            all_faces = union(obstacle_faces, free_faces)
            
            cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces);
            feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces));
            merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
            (valid_cover, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
            (valid_cover_merged, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover);
            # println("Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")

            if valid_cover && valid_cover_merged
                write(f,"Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")
                percent = round( ((length(cover)-length(merged_cover)) / length(cover)) * 100, digits = 2)
                write(f," | Percent Decrease: $percent\n")
                push!(percent_vec, percent)
            else
                fails += 1
                push!(tuples, (seed, num_obs))
                if !valid_cover_merged
                    write(f,"Merged cover invalid.\n")
                else
                    write(f,"Original cover invalid.\n")
                end
            end
        end
    end

    write(f, "\nTotal Fails: $fails.\nTuples: $tuples")
    flush(f)
    close(f)

    return percent_vec
end

#------------------------------------------------------------------------------------------------------

# method <- "merged" for the compact biclique cover, "full" for the original biclique cover, "bigM" for big-M constraints
function solve_time_stats(seed_range, num_obs_range; method="merged", file_name="Solve Time Stats.txt")
    # file_name = "Solve Time Stats.txt" #Seed Range = $seed_range, Num Obs Range = $num_obs_range.txt"
    # f = open(file_name, "a")
    # write(f, "Method = $method, Seed Range = $seed_range, Num Obs Range = $num_obs_range\n")
    # flush(f)
    # close(f)

    # Set parameters
    N = 30  # number of steps
    f1 = [0.0, 0.1, 0.0]  # initial footstep pose 1
    f2 = [0.0, 0.0, 0.0]  # initial footstep pose 2
    goal = [1, 1, 0]  # goal pose
    Q_g = 10*Matrix{Float64}(I, 3, 3)  # weight between final footstep and goal pose
    Q_r = [1 0 0; 0 1 0; 0 0 0] # Matrix{Float64}(I, 3, 3)  # weight between footsteps
    q_t = -.05  # weight for trimming unused steps
    # Optional named arguments
    # d1 = 0.2 <- radius of reference foot circle
    # d2 = 0.2 <- radius of moving foot circle
    # p1 = [0, 0.07] <- center of reference foot circle
    # p2 = [0, -0.27] <- center of moving foot circle
    # delta_x_y_max = 0.10  # max stride norm in space
    # delta_θ_max = pi/4  # max difference in θ
    # L = 5  # number of pieces in p.w.l approx. of sine/cosine
    
    f = open(file_name, "a")
    # write(f, "Seed\tNum_obs\tTime\n")
    write(f, "Seed\tNum_obs\tSolve_time\tRel_gap\tSimplex_iterations\tNodes_explored\tNum_vertices\tMerged_size\tFull_size\tNum_free_faces\tNum_free_face_inequalities\n")
    times = zeros(length(seed_range) * length(num_obs_range))
    i = 1
    for seed = seed_range
        for num_obs = num_obs_range
            # f = open(file_name, "a")
            # write(f, "Seed = $seed, Num_Obs = $num_obs\n")
            # flush(f)
            # close(f)
            # println("On test Seed = $seed, Num_Obs = $num_obs")

            # Create obstacles
            obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)

            x, y, θ, t, stats = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method=method)
            term_status, r_solve_time, rel_gap, simplex_iters, r_node_count, num_vertices, merged_size, full_size, num_free_faces, num_free_face_ineq, method_used = stats

            if method_used == method
                times[i] = r_solve_time # time
            else
                times[i] = -1
            end
            write(f, "$seed\t$num_obs\t$(times[i])\t$rel_gap\t$simplex_iters\t$r_node_count\t$num_vertices\t$merged_size\t$full_size\t$num_free_faces\t$num_free_face_ineq\n")
            flush(f)

            i += 1

            # # Trim excess steps
            # num_to_trim = length(filter(tj -> tj > 0.5, t[3:end]))
            # if num_to_trim % 2 == 0
            #     x = vcat(x[1:2], x[num_to_trim + 3 : end]);
            #     y = vcat(y[1:2], y[num_to_trim + 3 : end]);
            #     θ = vcat(θ[1:2], θ[num_to_trim + 3 : end]);
            # else
            #     x = vcat(x[1], x[num_to_trim + 3 : end]);
            #     y = vcat(y[1], y[num_to_trim + 3 : end]);
            #     θ = vcat(θ[1], θ[num_to_trim + 3 : end]);
            # end
            # plot_steps(obstacles, x, y, θ)
            # plot_circles(x, y, θ, p1=p1, p2=p2)

        end
    end

    # write(f, "\nEND")
    # flush(f)
    close(f)

    return times

end

seed_range = 21:30 # done 1:20
num_obs_range = 1:4
methods = ("merged", "full", "bigM")
times = Dict()
# method = "merged"
for method in methods
    file_name = "Solve Time Stats Method $method.txt"
    f = open(file_name, "a")   # write or append appropriately
    write(f, "Seed Range = $seed_range, Num Obs Range = $num_obs_range\nMethod = $method\n")
    flush(f)
    close(f)
    time = solve_time_stats(seed_range, num_obs_range, method=method, file_name=file_name)
    times[method] = time
end

file_name = "Solve Time Stats Comparison.txt" #Seed Range = $seed_range, Num Obs Range = $num_obs_range.txt"
f = open(file_name, "a")    # write or append appropriately
times_m_b = times["bigM"] .- times["merged"]
times_m_f = times["full"] .- times["merged"]
times_f_b = times["bigM"] .- times["full"]
write(f, "Seed Range = $seed_range, Num Obs Range = $num_obs_range\n")
write(f, "Merged_vs_BigM\tMerged_vs_Full\tFull_vs_BigM\n")
for i = 1:length(times_m_b)
    write(f, "$(times_m_b[i])\t$(times_m_f[i])\t$(times_f_b[i])\n")
end
flush(f)
close(f)


# Debugging Functions ---------------------------------------------------------------

# Check if free_faces has any duplicates---------------------------------------------
function face_duplicates(faces)
    dup = 0
    dup_ind = []
    col_faces = collect(faces)
    for i in 1:length(col_faces)
        for j = (i+1):length(col_faces)
            sort_i = sort(col_faces[i])
            sort_j = sort(col_faces[j])
            if sort_i == sort_j
                dup += 1
                push!(dup_ind, (i,j))
            end
        end
    end
    return dup, dup_ind
end

# Check if points has any duplicates---------------------------------------------
function point_duplicates(points)
    dup = 0
    dup_ind = []
    tol = Rational(1.e-13)
    for i in 1:(length(points))
        for j in (i+1):(length(points))
            #if (points[i].first == points[j].first) && (points[i].second == points[j].second)
            if abs(points[i].first - points[j].first) < tol && abs(points[i].second - points[j].second) < tol
                dup += 1
                push!(dup_ind, (i,j))
            end
        end
    end
    # for (i,j) in dup_ind
    #     println("Points $i and $j: \n ($(points[i])) \n ($(points[j]))")
    # end
    return dup, dup_ind
end

# Check if halfspaces in mapped have any duplicated points within them ---------
function points_in_hs_duplicates(mapped)
    dup = 0
    dup_ind = []
    for (l,hs) in enumerate(mapped)
        for i in 1:(length(hs))
            for j in (i+1):(length(hs))
                if (hs[i].first == hs[j].first) && (hs[i].second == hs[j].second)
                    dup += 1
                    push!(dup_ind, (l,i,j))
                end
            end
        end
    end
    return dup, dup_ind
end

# Check if any halfspaces in mapped are duplicates-----------------------------
function hs_duplicates(mapped)
    dup = 0
    dup_ind = []
    for i = 1:length(mapped)
        for j = (i+1):length(mapped)
            if mapped[i] == mapped[j]
                dup += 1
                push!(dup_ind, (i,j))
            end
        end
    end
    return dup, dup_ind
end

# Check to see if points are in order in mapped---------------------------------
function points_in_mapped_ordered(obstacles, mapped)
    for hs in mapped
        plot(legend=false,xlims=(-0.05,1.05),ylims=(-0.05,1.05));
        ClutteredEnvPathOpt.plot_field(obstacles)
        ClutteredEnvPathOpt.plot_lines(obstacles)
        for point in hs
            display(scatter!([point.first], [point.second]))
        end
    end
end

# Plot points individually----------------------------------------------------------
function plot_intersections_indiv(intersections)
    x = map(point -> point[1], intersections)
    y = map(point -> point[2], intersections)
    Plots.plot(legend=false,xlims=(-0.05,1.05),ylims=(-0.05,1.05))
    for i = 1:length(intersections)
        # Plots.plot(legend=false,xlims=(-0.05,1.05),ylims=(-0.05,1.05))
        display(Plots.scatter!([x[i]],[y[i]],title="Point $i:\n$(x[i])\n$(y[i])"))
    end
end

# Check if neighbors are correct----------------------------------------------------
function plot_neighbors(obstacles, neighbors, points)
    for (i,point) in enumerate(points)
        plot();
        ClutteredEnvPathOpt.plot_field(obstacles)
        ClutteredEnvPathOpt.plot_lines(obstacles)
        display(scatter!([point.first],[point.second],title="Point $i"))
        x_i_neighbors = map(node -> points[node].first, neighbors[i])
        y_i_neighbors = map(node -> points[node].second, neighbors[i])
        display(scatter!(x_i_neighbors,y_i_neighbors, title="Point $i. Neighbors $(neighbors[i])"))
        #println("Point $i has neighbors $(neighbors[i]).")
    end
end

function list_neighbors(neighbors, points)
    for i in 1:length(points)
        println("Point $i has neighbors $(neighbors[i]).")
    end
end

function duplicates(n)
    dup = 0
    for i in 1:length(n)
        for j in (i+1):length(n)
            if n[i] == n[j]
                dup += 1
            end
        end
    end
    return dup != 0
end

function suspect_neighbors(neighbors)
    suspect = []
    for n in 1:length(neighbors)
        if length(neighbors[n]) > 4 || duplicates(neighbors[n])
        #if duplicates(neighbors[n])
            push!(suspect, n)
        end
    end
    for point in suspect
        println("Point $point has neighbors $(neighbors[point])")
    end
    return suspect
end

# # Face finder ----------------------------------------------------------------------------
# t_seed = 7
# t_num_obs = 4
# obstacles = ClutteredEnvPathOpt.gen_field(t_num_obs, seed = t_seed) # optional seed
# construct_graph_debug(obstacles)