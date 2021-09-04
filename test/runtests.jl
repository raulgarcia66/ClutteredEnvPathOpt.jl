#ENV["GUROBI_HOME"] = "/Library/gurobi910/mac64"

using ClutteredEnvPathOpt
using LinearAlgebra
using Test
using Pipe
using Plots
using JuMP, Gurobi


# @testset "ClutteredEnvPathOpt.jl" begin
    # Create obstacles
    num_obs = 1;
    seed = 5;
    obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)
    plot();
    ClutteredEnvPathOpt.plot_field(obstacles)
    display(plot!())
    
    # Set parameters
    N = 30  # number of steps
    f1 = [0, 0.07, 0]  # initial footstep pose 1
    f2 = [0.0, 0, 0]  # initial footstep pose 2
    goal = [1, 1, 0]  # goal pose
    Q_g = 10*Matrix(I, 3, 3)  # weight between final footstep and goal pose
    Q_r = Matrix(I, 3, 3)  # weight between footsteps
    q_t = -.05  # weight for trimming unused steps
    L = 10  # number of pieces in p.w.l approx. of sine/cosine 
    delta_f_max = 1  # max stride norm
    # TODO: Assure second footstep is within reachability given first
    
    # Optional named arguments
    d1 = 0.20
    d2 = 0.20
    p1 = [0. 0.05]
    p2 = [0, -0.25]

    # Compute optimal path
    x, y, θ, t, z = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, L, delta_f_max)
    x2, y2, θ2, t2, z2 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, L, delta_f_max)
    # x, y, θ, t = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, L, delta_f_max, d1=d1, d2=d2, p1=p1, p2=p2)

    # Trim excess steps
    num_to_trim = length(filter(tj -> tj > 0.5, t[3:end]))
    x = vcat(x[1:2], x[num_to_trim + 3 : end]);
    y = vcat(y[1:2], y[num_to_trim + 3 : end]);
    θ = vcat(θ[1:2], θ[num_to_trim + 3 : end]);
    num_to_trim2 = length(filter(tj -> tj > 0.5, t2[3:end]))
    x2 = vcat(x2[1:2], x2[num_to_trim + 3 : end]);
    y2 = vcat(y2[1:2], y2[num_to_trim + 3 : end]);
    θ2 = vcat(θ2[1:2], θ2[num_to_trim + 3 : end]);

    # Plot footstep plan
    plot_steps(obstacles, x, y, θ)
    plot_steps(obstacles, x2, y2, θ2)
    png("Path Seed $seed Num Obs $num_obs")


    # Plot intersections of circles
    plot_circles(d1, d2, p1, p2, x2, y2, θ2)

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
num_obs = 3;
seed = 5;
obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs, "Obstacles Seed $seed Num Obs $num_obs", seed = seed)
skeleton = LabeledGraph(g)
all_faces = union(obstacle_faces, free_faces)

# Plot obstacles
plot();
ClutteredEnvPathOpt.plot_field(obstacles)
ClutteredEnvPathOpt.plot_lines(obstacles)
ClutteredEnvPathOpt.plot_borders()
ClutteredEnvPathOpt.plot_intersections(obstacles)
png("Obstacles with Lines Seed $seed Num Obs $num_obs")
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

feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces))
ClutteredEnvPathOpt.plot_edges(feg, points, obstacles, plot_name = "Finite Element Graph")
png("FEG Seed $seed Num Obs $num_obs")

feg_S = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
ClutteredEnvPathOpt.plot_edges(feg_S, points, obstacles, plot_name = "Finite Element Graph of S")
png("FEG of S Seed $seed Num Obs $num_obs")

feg_obs = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(obstacle_faces))
ClutteredEnvPathOpt.plot_edges(feg_obs, points, obstacles, plot_name = "Finite Element Graph Obstacle Faces")
png("FEG Obstacles Seed $seed Num Obs $num_obs")

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
ClutteredEnvPathOpt.plot_edges(missing_edges_lg, points, obstacles, plot_name = "Missing Edges")
png("Missing Edges Seed $seed Num Obs $num_obs Free")

file_name = "Biclique Cover Debug Output Seed $seed Num Obs $num_obs Free Unedited.txt"
f = open(file_name, "w")
write(f, "Plot 1\n")
flush(f)
close(f)
cover4 = ClutteredEnvPathOpt.find_biclique_cover_debug(skeleton, free_faces, points, obstacles, file_name)
# f = open("Biclique Cover Debug Output.txt","w")
# write(f, "Plot 1\n")
# flush(f)
# close(f)
# cover2 = ClutteredEnvPathOpt.find_biclique_cover_debug(skeleton, free_faces, points, obstacles)

# DON'T NEED TO ALL FACES COVER; JUST GIVING THE FREE FACES IS VALID
# All Faces Cover
# cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, all_faces)
# (valid_cover, size_diff, missing_edges, missing_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces)), cover)
# ClutteredEnvPathOpt.plot_edges(missing_edges_lg, points, obstacles, plot_name = "Missing Edges")
# png("Missing Edges Seed $seed Num Obs $num_obs All")

# file_name = "Biclique Cover Debug Output Seed $seed Num Obs $num_obs All Unedited.txt"
# f = open(file_name, "w")
# write(f, "Plot 1\n")
# flush(f)
# close(f)
# cover4 = ClutteredEnvPathOpt.find_biclique_cover_debug(skeleton, all_faces, points, obstacles, file_name)
# f = open("Biclique Cover Debug Output.txt","w")
# write(f, "Plot 1\n")
# flush(f)
# close(f)
# cover2 = ClutteredEnvPathOpt.find_biclique_cover_debug(skeleton, all_faces, points, obstacles)

# Biclique Cover Validity Tests
seed_range = 1:100
num_obs_range = 1:4
fails, tuples = biclique_cover_validity_tests(seed_range, num_obs_range, faces = "all")
# tuples = [(26,4), (100,4)] for seed_range = 1:100 and num_obs_range = 1:4
fails, tuples = biclique_cover_validity_tests(seed_range, num_obs_range, faces = "free")
# tuples = [(71,4)] for seed_range = 1:100 and num_obs_range = 1:4


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
                (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
            else
                cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces);
                feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces));
                (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
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