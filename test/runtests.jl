#ENV["GUROBI_HOME"] = "/Library/gurobi910/mac64"

using ClutteredEnvPathOpt
using LinearAlgebra
using Test
using Pipe
using Plots
using JuMP, Gurobi


# @testset "ClutteredEnvPathOpt.jl" begin
    obstacles = ClutteredEnvPathOpt.gen_field(2) # optional seed argument, default is 11

    x, y, theta = solve_deits(
        obstacles,
        20, # number of steps
        [0, 0, 0], # initial footstep 1
        [0.05, 0.05, 0], # initial footstep 2
        [1, 1, 0], # goal position to reach for footstep N
        Matrix(I, 3, 3),
        Matrix(I, 3, 3),
        -.0001, # q_t, must be negative
        10, # number of pieces in p.w.l approx. of sine/cosine
        1 #
    )
    plot();
    ClutteredEnvPathOpt.plot_field(obstacles)
    #ClutteredEnvPathOpt.plot_lines(obstacles)
    #ClutteredEnvPathOpt.plot_intersections(obstacles)
    display(plot!(xlims=(-0.05,1.05), ylims= (-0.05,1.05)));

    scatter!(x[1:2:end], y[1:2:end], color="red", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:2:length(x)]))
    scatter!(x[2:2:end], y[2:2:end], color="blue", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 2:2:length(x)]))
    quiver!(x, y, quiver=(0.1 * cos.(theta), 0.1 * sin.(theta)))

#     @test 1 == 1
# end

# @testset "Optimal biclique comparison" begin
    obstacles, points, g, faces = ClutteredEnvPathOpt.plot_new(2,"Biclique Comparison Test") # optional seed argument, default is 11
    skeleton = LabeledGraph(g)

    tree = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton, faces)
    (cover_vec, graphviz) = ClutteredEnvPathOpt.tree2digraph(tree)
    io = open("graphviz.txt", "a")
    println(io, graphviz)
    close(io)
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
#     lower = Int(ceil(log2(length(faces))))
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


# Unofficial Tests ---------------------------------------------------------------------------

# There is a paradigm for constructing proper tests
@testset "Biclique Edge Visualization" begin
    num_obs = 4;
    seed = 9;
    obstacles, points, g, free_faces = ClutteredEnvPathOpt.plot_new(num_obs,"Obstacles Seed $seed Num Obs $num_obs",seed) # optional seed argument, default is 11
    png("All Faces Seed $seed Num Obs $num_obs")
    skeleton = LabeledGraph(g)

    plot();
    ClutteredEnvPathOpt.plot_field(obstacles)
    ClutteredEnvPathOpt.plot_lines(obstacles)
    ClutteredEnvPathOpt.plot_intersections(obstacles)
    display(plot!(legend=false))

    plot();
    plot_free_faces(free_faces)
    png("Free Faces Seed $seed Num Obs $num_obs")
    dup, dup_ind = face_duplicates(free_faces)

    plot_faces_indiv_with_points(free_faces, points)

    plot_edges(skeleton, points, obstacles)
    png("Skeleton Seed $seed Num Obs $num_obs")
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
    plot_edges(feg, points, obstacles)

    obstacles = ClutteredEnvPathOpt.gen_field(num_obs,seed) # optional seed
    points, mapped = ClutteredEnvPathOpt.find_intersections(obstacles)
    neighbors = ClutteredEnvPathOpt._find_neighbors(points, mapped)

    plot_intersections_indiv(points)
    dup, dup_ind = point_duplicates(points)
    dup, dup_ind = points_in_hs_duplicates(mapped)
    dup, dup_ind = hs_duplicates(mapped)
    points_in_mapped_ordered(obstacles, mapped)
    plot_neighbors(obstacles, neighbors, points)
    list_neighbors(neighbors, points)
    suspect = suspect_neighbors(neighbors)

    # construct_graph_debug(obstacles)

    #cover = ClutteredEnvPathOpt._find_biclique_cover_visual(skeleton, free_faces, points)
    cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces)
    # cover_vec = collect(cover)
    ClutteredEnvPathOpt._is_valid_biclique_cover_diff(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces)), cover)
end


# Debugging Functions ------------------------------------------------------------
import Polyhedra
import Statistics
import LightGraphs
import GLPK

# Plot free faces ----------------------------------------------------------------
function plot_free_faces(free_faces)
    for (j,face) in enumerate(free_faces)
        #println("$face")
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        x_locations = map(i -> points[i].first, face)
        y_locations = map(i -> points[i].second, face)
        #println("$x_locations")

        avg_x = Statistics.mean(x_locations)
        avg_y = Statistics.mean(y_locations)
        #println("$avg_x , $avg_y")
        polygon = Polyhedra.polyhedron(
            v,
            Polyhedra.DefaultLibrary{Float64}(Gurobi.Optimizer)
        )

        plot!(polygon)
        plot!([avg_x], [avg_y], series_annotations=([Plots.text("$j", :center, 8, "courier")]))
        display(plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))

        # display(plot!(polygon, title="Free Face # $j", xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
    end
    # display(plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
end

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

# Plot faces indiviually with their points to check if clockwise-ordered
function plot_faces_indiv_with_points(faces, points)
    for (j,face) in enumerate(faces)
        plot()
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        x_locations = map(i -> points[i].first, face)
        y_locations = map(i -> points[i].second, face)

        # avg_x = Statistics.mean(x_locations)
        # avg_y = Statistics.mean(y_locations)
        # println("$avg_x , $avg_y")
        polygon = Polyhedra.polyhedron(
            v,
            Polyhedra.DefaultLibrary{Float64}(Gurobi.Optimizer)
        )

        plot!(polygon)
        plot!(x_locations, y_locations, series_annotations=([Plots.text(string(x), :center, 8, "courier") for x in face]))
        display(plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05),title="Face $j: $face"))
        # display(plot!(title="Face $j: $face"))

    end
end

# Check if points has any duplicates---------------------------------------------
function point_duplicates(points)
    dup = 0
    dup_ind = []
    tol = 1.e-13
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

# Check if halfspaces in mapped have any duplicated points within the halfspace---
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

# Plot points individually
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

# Plot edges of a lightgraph
function plot_edges(lg, points, obstacles)
    Plots.plot(xlims=(-0.05,1.05),ylims=(-0.05,1.05),legend=false)
    ClutteredEnvPathOpt.plot_intersections(obstacles)
    for edge in LightGraphs.edges(lg.graph)
        Plots.plot!([points[edge.src].first, points[edge.dst].first], [points[edge.src].second, points[edge.dst].second])
        display(Plots.plot!(title="Edge ($(edge.src), $(edge.dst))"))
    end
    #display(Plots.plot!(title="Edges"))
end

# Face finder
t_seed = 7
t_num_obs = 4
obstacles = ClutteredEnvPathOpt.gen_field(t_num_obs,t_seed) # optional seed
construct_graph_debug(obstacles)

function construct_graph_debug(obs)
    points, mapped = ClutteredEnvPathOpt.find_intersections(obs)

    # Create map from point to neighbors (counter clockwise ordered by angle against horizontal)
    neighbors = ClutteredEnvPathOpt._find_neighbors(points, mapped)
    graph = ClutteredEnvPathOpt._gen_graph(neighbors)

    angles = fill(Rational(Inf), length(points), length(points))
    for i in 1:length(points)
        for j in 1:length(points)
            angles[i, j] = Rational(atan(points[j].second - points[i].second, points[j].first - points[i].first))
            if angles[i, j] < 0//1
                angles[i, j] += 2//1 * Rational(3.141592653589793) # * pi
            end
        end
    end

    function greatest_angle_neighbor(source, last_angle, visited)
        #tol = Rational(1.e-15)
        unvisited = filter(vertex -> !(vertex in visited), neighbors[source])
        display(println("Unvisisted nodes are $unvisited"))
        neighbor_angles = map(neighbor -> angles[source, neighbor], unvisited)
        for i in 1:length(neighbor_angles)
            if neighbor_angles[i] >= Rational(3.141592653589793)
                neighbor_angles[i] -= Rational(2*3.141592653589793)
            end
        end
        zipped = @pipe zip(unvisited, neighbor_angles) |> collect(_)
        display(map(zip -> println("Unvisited nodes and angles are ($(zip[1]),($(Float64(zip[2])))"), zipped))
        # zipped # zipped is an array of tuples of the form (unvisited node, angle it makes with source)

        # first condition is to assure our face is the most compact one
        # second condition prevents selecting a node on the same halfspace as source
        #greater = filter(tup -> ((tup[2] > last_angle) && !(abs(tup[2] - last_angle) < tol)), zipped)
        # greater = filter(tup -> (tup[2] > last_angle), zipped)

        if last_angle >= 0 && last_angle < Rational(3.141592653589793)
            greater = filter(tup -> (tup[2] > last_angle) && (tup[2] < last_angle + Rational(3.141592653589793)), zipped)
        else
            greater = filter(tup -> (tup[2] > last_angle - Rational(2*3.141592653589793)) && (tup[2] < last_angle - Rational(3.141592653589793)), zipped)
        end

        if isempty(greater)
            return (-1, -1) # bad state
        else
            sort!(greater, by=(tup -> tup[2]))
            display(map(zip -> println("greater: nodes and angles are ($(zip[1]),($(Float64(zip[2])))"), greater))
            destination = greater[end][1] # greater[1][1]
        end

        return (destination, angles[source, destination])
    end

    faces = []

    for i in 1:length(points)
    # for i in 7:10
        start = i
        println("\nStart is $start -----------------------------------------------------------")

        for j in 1:length(neighbors[start])
            current = neighbors[start][j]
            last_angle = angles[start, current]
            println("\nCurrent is $current. last_angle is $(Float64(last_angle))")

            face = [start, current]

            while start != current
                # current # ?
                println("Face length: $(length(face)). Start in neighbors[current] $(start in neighbors[current])")
                if (length(face) > 2 && (start in neighbors[current]))
                    push!(faces, copy(face))
                    println("$start is in the neighbors of $current")
                    println("Face pushed is $face\n########")
                    break
                end

                current, last_angle = greatest_angle_neighbor(current, last_angle, face)
                println("Update: Current is now $current. last_angle is now $(Float64(last_angle))")

                if current == -1
                    break
                end

                push!(face, current)
                println("Update: Face is now $face")
            end
        end
    end

    # Cleaning
    face_sets = []

    # Get the nodes of the extreme points of the obstacles
    # obstacle_faces = map(
    #     obstacle -> Set(
    #         map(
    #             point -> findfirst(p -> (isapprox(p.first, point[1]) && isapprox(p.second, point[2])), points),
    #             obstacle.vrep.points.points
    #         )
    #     ),
    #     obs
    # )

    # Removes redundant faces and faces that comprise an obstacle
    unique_faces = filter(face -> begin
        face_set = Set(face)

        face_set_v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        # face_set_poly = Polyhedra.polyhedron(face_set_v, Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer))
        face_set_poly = Polyhedra.polyhedron(face_set_v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
        face_set_overlaps_obs_faces = false
        included = false

        tol =  Rational(1.e-10)
        for n in 1:length(obs) # need to loop in this manner for 'continue' to work properly
            indicator = 0
            intersect_poly = Polyhedra.intersect(face_set_poly, obs[n])
            Polyhedra.npoints(intersect_poly) # computes the extreme points so that vrep can be used

            if typeof(intersect_poly.vrep) === nothing
                continue
            end
            ext_p_intersect_poly = sort(collect(intersect_poly.vrep.points.points))
            ext_p_face_set_poly = sort(collect(face_set_poly.vrep.points.points))
            if length(ext_p_intersect_poly) != length(ext_p_face_set_poly)
                continue
            end
            for i in 1:length(ext_p_face_set_poly)
                for j in 1:length(ext_p_face_set_poly[i])
                    if abs(ext_p_intersect_poly[i][j] - ext_p_face_set_poly[i][j]) > tol
                        indicator += 1
                        break
                    end
                end
                if indicator != 0
                    break
                end
            end
            if indicator == 0
                face_set_overlaps_obs_faces = true
                break
            end
        end

        if !(face_set in face_sets) && !(face_set_overlaps_obs_faces)
            push!(face_sets, face_set)
            included = true
        end

        return included
    end, faces)

    # Plot all faces to see if all are found
    Plots.plot()
    for (j,face) in enumerate(faces)
        #Plots.plot()
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        x_locations = map(i -> points[i].first, face)
        y_locations = map(i -> points[i].second, face)
        #println("$x_locations")

        avg_x = Statistics.mean(x_locations)
        avg_y = Statistics.mean(y_locations)
        #println("$avg_x , $avg_y")
        polygon = Polyhedra.polyhedron(
            v,
            # Polyhedra.DefaultLibrary{Float64}(Gurobi.Optimizer)
            Polyhedra.DefaultLibrary{Rational{Int64}}(Gurobi.Optimizer)
        )

        Plots.plot!(polygon, title="$j-th Face: $face")
        Plots.plot!([avg_x], [avg_y], series_annotations=([Plots.text("$j", :center, 8, "courier")]))
        display(Plots.plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))

    end
    # display(Plots.plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))

    return (obs, points, graph, @pipe map(face -> reverse(face), unique_faces) |> Set(_))
    # return (obs, points, graph, @pipe map(face -> reverse(face), faces) |> Set(_))
end