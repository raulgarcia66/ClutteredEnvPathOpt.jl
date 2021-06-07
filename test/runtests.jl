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
    obstacles, points, g, faces = ClutteredEnvPathOpt.plot_new(2,"B") # optional seed argument, default is 11
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
    obstacles, points, g, free_faces = ClutteredEnvPathOpt.plot_new(2,"Biclique Obstacles Seed 7",7) # optional seed arument, default is 11
    skeleton = LabeledGraph(g)
    
    plot();
    plot_free_faces(free_faces) # Need to run this function below first

    cover = ClutteredEnvPathOpt._wrapper_find_biclique_cover(skeleton, free_faces, points)
    # cover_vec = collect(cover)
    ClutteredEnvPathOpt._is_valid_biclique_cover(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces)), cover)
end


# Plot free spaces and obstacles -----------------------------------------------
import Polyhedra
import Statistics

obstacles, points, g, free_faces = ClutteredEnvPathOpt.plot_new(2,"Biclique Obstacles #") # optional seed

plot()
intersections = ClutteredEnvPathOpt.find_intersections(obstacles)[1]
x = map(point -> point[1], intersections)
y = map(point -> point[2], intersections)
scatter!(x,y, series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:length(x)]))

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
        #display(plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))

        #display(plot(polygon, title="Free Face # $j", xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
    end
    display(plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
end
#plot_free_faces(free_faces)

# Check if free_faces has any duplicates---------------------------------------------
dup = 0
dup_ind = (-1,-1)
col_free_faces = collect(free_faces)
for i in 1:length(col_free_faces)
    for j = (i+1):length(col_free_faces)
        sort_i = sort(col_free_faces[i])
        sort_j = sort(col_free_faces[j])
        if sort_i == sort_j
            dup += 1
            dup_ind = (i,j)
        end
    end
end
dup
dup_ind


# Create obstacles and related data ---------------------------------------------------
plot();
obstacles = ClutteredEnvPathOpt.gen_field(2)
points, mapped = ClutteredEnvPathOpt.find_intersections(obstacles)
ClutteredEnvPathOpt.plot_field(obstacles)
#ClutteredEnvPathOpt.plot_lines(obstacles)
# To plot intersecions
# intersections = ClutteredEnvPathOpt.find_intersections(obstacles)[1]
# x = map(point -> point[1], intersections)
# y = map(point -> point[2], intersections)
# scatter!(x,y, series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:length(x)]))
display(plot!())

# Check if points has any duplicates---------------------------------------------
dup = 0
for i in 1:(length(points)-1)
    if points[i].first == points[i+1].first
        if points[i].second == points[i+1].second 
            duplicates += 1
        end
    elseif points[i].first == points[i+1].second 
        if points[i].second == points[i+1].first
            duplicates += 1
        end
    end
end
dup

# Check to see if the points are in order-------------------------------------------
for hs in mapped
    for point in hs
        display(scatter!([point.first], [point.second]))
    end
end

# Check if neighbors are correct----------------------------------------------------
neighbors = ClutteredEnvPathOpt._find_neighbors(points, mapped)
plot()
ClutteredEnvPathOpt.plot_field(obs)
display(plot!())
for (i,point) in enumerate(points)
    display(scatter!([point.first],[point.second],title="Point $i"))
    x_i_neighbors = map(node -> points[node].first, neighbors[i])
    y_i_neighbors = map(node -> points[node].second, neighbors[i])
    display(scatter!(x_i_neighbors,y_i_neighbors))
end
suspect = []
for n in 1:length(neighbors)
    if length(neighbors[n]) < 4
        push!(suspect, n)
    end
end
for point in suspect
    println("Point $point has neighbors $(neighbors[point])")
end

# Test overlap detection procedure-----------------------------------------------------------------
# Need faces, obs, points
counter = 0
for (n,face) in enumerate(faces)
    # faces = collect(free_faces)
    #face = faces[9]
    face_set = Set(face)
    #included = face_set in face_sets

    face_set_v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
    face_set_poly = Polyhedra.polyhedron(face_set_v, Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer))
    face_set_overlaps_obs_faces = false
    included = false

    tol =  1.e-6
    #for (i,obs_face) in enumerate(obs)
    #for obs_face in obs
    for l in 1:length(obs) # need to loop like this for "continue" to work properly
        indicator = 0
        intersect_poly = Polyhedra.intersect(face_set_poly,obs[l]) #obs_face)
        Polyhedra.npoints(intersect_poly)
        # try
        #     intersect_poly.vrep.points.points 
        # catch e
        #     println("No v points. Obstacle 2")
        # end

        if typeof(intersect_poly.vrep) == Nothing
            println("Face $n: Intersection is empty with obstacle $l")
            continue
        end
        ext_p_intersect_poly = sort(collect(intersect_poly.vrep.points.points))
        ext_p_face_set_poly = sort(collect(face_set_poly.vrep.points.points))
        if length(ext_p_intersect_poly) != length(ext_p_face_set_poly)
            println("Different number of extreme points")
            continue
        end
        for i in 1:length(ext_p_face_set_poly)
            for j in 1:length(ext_p_face_set_poly[i])
                #println("$(ext_p_intersect_poly[i][j]) - $(ext_p_face_set_poly[i][j]) = $(ext_p_intersect_poly[i][j] - ext_p_face_set_poly[i][j])")
                if abs(ext_p_intersect_poly[i][j] - ext_p_face_set_poly[i][j]) > tol
                    println("Tolerance breach")
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
            counter += 1
            println("$face_set_overlaps_obs_faces")
            break
        end
    end
end