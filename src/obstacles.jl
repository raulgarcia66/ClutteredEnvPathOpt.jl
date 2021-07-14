import LightGraphs
import Polyhedra
import GLPK
import Plots
import Random

import Statistics
import Gurobi

using Pipe

"""
    gen_obstacle(max_side_len)

Creates a triangle at a random location in the unit square with given max side
length.
"""
function gen_obstacle(max_side_len::Real)
    point = rand(2)

    v = Polyhedra.convexhull(point, map(x -> x > 1 ? 1 : x, point + max_side_len * rand(2)), map(x -> x > 1 ? 1 : x, point + max_side_len * rand(2)))
    
    # Non-random; function will create this obstacle every time it's called
    # v = Polyhedra.convexhull([.25;.25], [.25;.5], [.5;.25])
    
    return Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer))
end

"""
    gen_obstacle_finite_points(points_per_dim, max_jumps)

Creates a random triangle from a subset of points in a grid, with a max distance between points.
"""
function gen_obstacle_finite_points(points_per_dim::Int, max_jumps::Int)
    nums = Rational.(LinRange(0//1,1//1,points_per_dim+2))
    nums = nums[2:end-1]
    x_index_1 = rand(1:points_per_dim)
    y_index_1 = rand(1:points_per_dim)
    point_1 = [ nums[x_index_1] ; nums[y_index_1] ]

    min_x_left = min(x_index_1 -1, max_jumps)
    min_x_right = min(points_per_dim - x_index_1, max_jumps)
    min_y_bottom = min(y_index_1 -1, max_jumps)
    min_y_top = min(points_per_dim - y_index_1, max_jumps)

    same_point_as_1 = true
    x_index_2 = -1
    y_index_2 = -1
    while same_point_as_1
        x_index_2 = rand( x_index_1 - min_x_left : x_index_1 + min_x_right )
        y_index_2 = rand( y_index_1 - min_y_bottom: y_index_1 + min_y_top )
        if x_index_2 != x_index_1 || y_index_2 != y_index_1
            same_point_as_1 = false
        end
    end
    point_2 = [ nums[x_index_2] ; nums[y_index_2] ]

    same_point_as_1_and_2 = true
    x_index_3 = -1
    y_index_3 = -1
    while same_point_as_1_and_2
        x_index_3 = rand( x_index_1 - min_x_left : x_index_1 + min_x_right )
        y_index_3 = rand( y_index_1 - min_y_bottom : y_index_1 + min_y_top )
        if (x_index_3 != x_index_1 || y_index_3 != y_index_1) && (x_index_3 != x_index_2 || y_index_3 != y_index_2)
            same_point_as_1_and_2 = false
        end
    end
    point_3 = [ nums[x_index_3] ; nums[y_index_3] ]

    v = Polyhedra.convexhull(point_1, point_2, point_3)
    # poly = Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer))
    poly = Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
    if length(poly.vrep.points.points) != 3
        return gen_obstacle_finite_points(points_per_dim, max_jumps)
    end

    # display(Plots.plot(poly, xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
    # return Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
    return poly
end

"""
    gen_field(num_obstacles)

Creates an obstacle course, an array of obstacles in the unit square. If
obstacles overlap the convex hull of their combined vertices will be taken
instead.
"""
function gen_field(num_obstacles::Int, seed::Int=11)
    Random.seed!(seed)
    # obstacles = map(_n -> gen_obstacle(0.25), 1:num_obstacles)
    obstacles = map(_n -> gen_obstacle_finite_points(10,3), 1:num_obstacles)

    return remove_overlaps(obstacles)
end

"""
    remove_overlaps(obstacles)

Given an obstacle course that may contain overlapping obstacles, check for
overlaps and replace offending obstacles with the convex hull of their combined
vertices.
"""
function remove_overlaps(obstacles)
    if !has_overlaps(obstacles)
        return obstacles
    end

    res = []
    fused = []

    for i = 1:length(obstacles)
        if i in fused
            continue
        end

        accumulator = obstacles[i]
        for j = 1:length(obstacles)
            if i != j && !(j in fused)
                other = obstacles[j]

                overlap = Polyhedra.intersect(accumulator, other)
                if !isempty(overlap)
                    accumulator = Polyhedra.convexhull(accumulator, other)
                    push!(fused, j)
                end
            end
        end

        push!(res, accumulator)
    end

    return remove_overlaps(res)
end

"""
    has_overlaps(obstacles)

Returns true if any obstacles in the set are overlapping.
"""
function has_overlaps(obstacles)::Bool
    for i = 1:length(obstacles)
        for j = 1:length(obstacles)
            if i != j
                if !isempty(Polyhedra.intersect(obstacles[i], obstacles[j]))
                    return true
                end
            end
        end
    end

    return false
end

"""
    find_intersections(obstacles)

Analyzes the halfspaces that make up each obstacle and finds where the
halfspaces intersect. Returns an array of all intersections as well as an
array mapping each halfspace to the points of intersection that lie on it in
order. Will also consider the lines the make up the border of the unit square.
"""
function find_intersections(obstacles)
    res = Set([])
    
    halfspaces = @pipe map(obstacle -> Polyhedra.hrep(obstacle).halfspaces, obstacles) |> Iterators.flatten(_) |> collect(_)
    unique!(halfspaces)

    boundaries = [
        Polyhedra.HalfSpace([-1//1, 0//1], 0//1),
        Polyhedra.HalfSpace([1//1, 0//1], 1//1),
        Polyhedra.HalfSpace([0//1, -1//1], 0//1),
        Polyhedra.HalfSpace([0//1, 1//1], 1//1)
    ]

    # unique_boundaries = []
    # for b in boundaries
    #     if !(b in halfspaces)
    #         push!(unique_boundaries, b)
    #     end
    # end
    # lines = [halfspaces; unique_boundaries]

    # Remove duplicate halfspaces
    filter!(b -> !(b in halfspaces), boundaries)
    lines = [halfspaces; boundaries]
    # for b in boundaries
    #     if (b in halfspaces)
    #         ind = findfirst(x -> x == b, boundaries)
    #         if ind !== nothing
    #             splice!(boundaries, ind)
    #         end
    #     end
    # end
    # lines = [halfspaces; boundaries]

    # Save indices of horizontal and vertical halfspaces
    vert_ind = []
    hor_ind = []
    for i = 1:length(lines)
        if abs(lines[i].a'[1]) == 1//1 && abs(lines[i].a'[2]) == 0//1
            push!(vert_ind, i)
        end
        if abs(lines[i].a'[1]) == 0//1 && abs(lines[i].a'[2]) == 1//1
            push!(hor_ind, i)
        end
    end

    points = fill([], length(lines))    # Association from lines to points
    for i = 1:length(lines)   # Does same as line above? Apparently not
        points[i] = []
    end

    function point_duplicates(point, res)
        tol = Rational(1.e-13)
        for index in 1:length(res)
            if abs(point.first - res[index].first) < tol && abs(point.second - res[index].second) < tol 
            # if point.first == res[index].first && point.second == res[index].second
                return true, index
            end
        end
        return false, -1
    end
    
    tol = Rational(1.e-13)
    x_vert = Dict([]) # keys are the indices
    y_hor = Dict([]) 
    for i = 1:length(lines)     
        for j = (i + 1):length(lines)
            A = vcat(lines[i].a', lines[j].a')
            b = [lines[i].β; lines[j].β]
            try
                x = A \ b

                if abs(x[1] - 0//1) < tol
                    x[1] = 0//1
                elseif abs(x[1] - 1//1) < tol
                    x[1] = 1//1
                end

                if abs(x[2] - 0//1) < tol
                    x[2] = 0//1
                elseif abs(x[2] - 1//1) < tol
                    x[2] = 1//1
                end

                # point = x[1] => x[2]

                if (x[1] >= 0//1 && x[1] <= 1//1 && x[2] >= 0//1 && x[2] <= 1//1)
                    # Assures points on horizontal or vertical halfspaces have the same y- or x- coordinates, resp.
                    if i in vert_ind
                        if !(i in keys(x_vert))
                            x_vert[i] = x[1]
                        else
                            x[1] = x_vert[i]
                        end
                    elseif i in hor_ind
                        if !(i in keys(y_hor))
                            y_hor[i] = x[2]
                        else
                            x[2] = y_hor[i]
                        end
                    end
                    if j in vert_ind
                        if !(j in keys(x_vert))
                            x_vert[j] = x[1]
                        else
                            x[1] = x_vert[j]
                        end
                    elseif j in hor_ind
                        if !(j in keys(y_hor))
                            y_hor[j] = x[2]
                        else
                            x[2] = y_hor[j]
                        end
                    end

                    point = x[1] => x[2]

                    dup, index = point_duplicates(point, collect(res))
                    if !(dup)
                        push!(res, point)
                        push!(points[i], point)
                        push!(points[j], point)
                    else
                        if !(collect(res)[index] in points[i])
                            push!(points[i], collect(res)[index])
                        end
                        if !(collect(res)[index] in points[j])
                            push!(points[j], collect(res)[index])
                        end
                    end
                end
            catch e
                # No intersections found
            end
        end
    end

    sorted_points = map(list -> sort(list), points)

    return (collect(res), sorted_points)
end

"""
    construct_graph(obs)

Given a list of obstacles with no overlaps finds the intersections of the lines
that make up each obstacles halfspace and constructs a graph with the
intersections as nodes and the halfspaces as edges. Returns a tuple containing
the obstacles, a list of the location of the nodes' coordinates, the graph, and
a set of every face of free space (where a face is a clockwise-ordered vector
of vertices).
"""
function construct_graph(obs)
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
        # display(println("Unvisisted nodes are $unvisited"))
        neighbor_angles = map(neighbor -> angles[source, neighbor], unvisited)
        # For sorting angles
        for i in 1:length(neighbor_angles)
            if neighbor_angles[i] >= Rational(3.141592653589793)
                neighbor_angles[i] -= Rational(2*3.141592653589793)
            end
        end
        zipped = @pipe zip(unvisited, neighbor_angles) |> collect(_)
        # display(map(zip -> println("Unvisited nodes and angles are ($(zip[1]),($(Float64(zip[2])))"), zipped))
        # zipped is an array of tuples of the form (unvisited node, angle it makes with source)

        # Two cases to assure we are constructing the most compact face (ie, no face with subfaces)
        if last_angle >= 0 && last_angle < Rational(3.141592653589793)
            greater = filter(tup -> (tup[2] > last_angle) && (tup[2] < last_angle + Rational(3.141592653589793)), zipped)
        else
            greater = filter(tup -> (tup[2] > last_angle - Rational(2*3.141592653589793)) && (tup[2] < last_angle - Rational(3.141592653589793)), zipped)
        end

        if isempty(greater)
            return (-1, -1) # bad state
        else
            sort!(greater, by=(tup -> tup[2]))
            #display(map(zip -> println("greater: nodes and angles are ($(zip[1]),($(Float64(zip[2])))"), greater))
            destination = greater[end][1] # greater[1][1]
        end

        return (destination, angles[source, destination])
    end

    faces = []

    for i in 1:length(points)
        start = i
        #println("\nStart is $start -----------------------------------------------------------")

        for j in 1:length(neighbors[start])
            current = neighbors[start][j]
            last_angle = angles[start, current]

            face = [start, current]

            while start != current

                if (length(face) > 2 && (start in neighbors[current]))
                    push!(faces, copy(face))
                    break
                end

                current, last_angle = greatest_angle_neighbor(current, last_angle, face)

                if current == -1
                    break
                end

                push!(face, current)
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
        # display(Plots.plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))

    end
    display(Plots.plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))

    return (obs, points, graph, @pipe map(face -> reverse(face), unique_faces) |> Set(_))
    # return (obs, points, graph, @pipe map(face -> reverse(face), faces) |> Set(_))
end

# function construct_graph(obs)
#     points, mapped = find_intersections(obs)
    
#     # Create map from point to neighbors (counter clockwise ordered by angle against horizontal)
#     neighbors = _find_neighbors(points, mapped)
#     graph = _gen_graph(neighbors)

#     # angles = map(point -> atan(point.first, point.second), points)
#     angles = fill(Rational(Inf), length(points), length(points))
#     for i in 1:length(points)
#         for j in 1:length(points)
#             angles[i, j] = Rational.(atan(points[j].second - points[i].second, points[j].first - points[i].first))
#             if angles[i, j] < 0//1
#                 angles[i, j] += 2//1 * Rational(3.141592653589793) # * pi              
#             end
#         end
#     end

#     function greatest_angle_neighbor(source, last_angle, visited)
#         tol = Rational(1.e-15)
#         unvisited = filter(vertex -> !(vertex in visited), neighbors[source])
#         neighbor_angles = map(neighbor -> angles[source, neighbor], unvisited)
#         zipped = @pipe zip(unvisited, neighbor_angles) |> collect(_)
#         zipped # zipped is an array of tuples of the form (unvisited node, angle it makes with source)

#         # greater = filter(tup -> ((tup[2] > last_angle) && !(isapprox(tup[2], last_angle))), zipped)
#         greater = filter(tup -> ((tup[2] > last_angle) && !(abs(tup[2] - last_angle) < tol)), zipped)


#         if isempty(greater)
#             return (-1, -1) # bad state
#         else
#             sort!(greater, by=(tup -> tup[2]))
#             destination = greater[1][1]
#         end

#         return (destination, angles[source, destination])
#     end

#     faces = []

#     for i in 1:length(points)
#         start = i

#         for j in 1:length(neighbors[start])
#             current = neighbors[start][j]
#             last_angle = angles[start, current]

#             face = [start, current]

#             while start != current
#                 current # ?
#                 if (length(face) > 2 && (start in neighbors[current]))
#                     push!(faces, copy(face))
#                     break;
#                 end

#                 current, last_angle = greatest_angle_neighbor(current, last_angle, face)

#                 if current == -1
#                     break
#                 end

#                 push!(face, current)
#             end
#         end
#     end

#     # Cleaning
#     face_sets = []

#     # Get the nodes of the extreme points of the obstacles
#     # obstacle_faces = map(
#     #     obstacle -> Set(
#     #         map(
#     #             point -> findfirst(p -> (isapprox(p.first, point[1]) && isapprox(p.second, point[2])), points),
#     #             obstacle.vrep.points.points
#     #         )
#     #     ),
#     #     obs
#     # )

#     # Removes redundant faces and faces that comprise an obstacle
#     unique_faces = filter(face -> begin
#         face_set = Set(face)

#         face_set_v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
#         # face_set_poly = Polyhedra.polyhedron(face_set_v, Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer))
#         face_set_poly = Polyhedra.polyhedron(face_set_v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
#         face_set_overlaps_obs_faces = false
#         included = false

#         tol =  Rational(1.e-10)
#         for n in 1:length(obs) # need to loop in this manner for 'continue' to work properly
#             indicator = 0
#             intersect_poly = Polyhedra.intersect(face_set_poly, obs[n])
#             Polyhedra.npoints(intersect_poly) # computes the extreme points so that vrep can be used

#             if typeof(intersect_poly.vrep) === nothing
#                 continue
#             end
#             ext_p_intersect_poly = sort(collect(intersect_poly.vrep.points.points))
#             ext_p_face_set_poly = sort(collect(face_set_poly.vrep.points.points))
#             if length(ext_p_intersect_poly) != length(ext_p_face_set_poly)
#                 continue
#             end
#             for i in 1:length(ext_p_face_set_poly)
#                 for j in 1:length(ext_p_face_set_poly[i])
#                     if abs(ext_p_intersect_poly[i][j] - ext_p_face_set_poly[i][j]) > tol
#                         indicator += 1
#                         break
#                     end
#                 end
#                 if indicator != 0
#                     break
#                 end
#             end
#             if indicator == 0
#                 face_set_overlaps_obs_faces = true
#                 break
#             end
#         end

#         if !(face_set in face_sets) && !(face_set_overlaps_obs_faces)
#             push!(face_sets, face_set)
#             included = true
#         end

#         return included
#     end, faces)

#     # Plot all faces to see if all are found
#     Plots.plot()
#     for (j,face) in enumerate(faces)
#         #Plots.plot()
#         v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
#         x_locations = map(i -> points[i].first, face)
#         y_locations = map(i -> points[i].second, face)
#         #println("$x_locations")
        
#         avg_x = Statistics.mean(x_locations)
#         avg_y = Statistics.mean(y_locations)
#         #println("$avg_x , $avg_y")
#         polygon = Polyhedra.polyhedron(
#             v,
#             # Polyhedra.DefaultLibrary{Float64}(Gurobi.Optimizer)
#             Polyhedra.DefaultLibrary{Rational{Int64}}(Gurobi.Optimizer)
#         )
    
#         Plots.plot!(polygon, title="$j-th Face: $face")
#         Plots.plot!([avg_x], [avg_y], series_annotations=([Plots.text("$j", :center, 8, "courier")]))
#         display(Plots.plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
    
#     end
#     # display(Plots.plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))

#     return (obs, points, graph, @pipe map(face -> reverse(face), unique_faces) |> Set(_))
# end

"""
    _find_neighbors(points, mapped)

Given a list of points and an association of halfspace lines to points that lie
on the line, return an array where res[i] <- the neighbors of node i
"""
function _find_neighbors(points, mapped)
    neighbors = Dict()

    for point in points
        neighbors[point] = []

        for colinear in mapped
            if in(point, colinear)
                index = findfirst(p -> p == point, colinear) # adjust for potential tol issue? Shouldn't arise

                if index - 1 > 0 && !(colinear[index - 1] in neighbors[point])
                    push!(neighbors[point], colinear[index - 1])
                end

                if index + 1 <= length(colinear) && !(colinear[index + 1] in neighbors[point])
                    push!(neighbors[point], colinear[index + 1])
                end
            end
        end
    end

    res = fill([], length(points))
    for i in 1:length(res)
        res[i] = map(point -> findfirst(p -> p == point, points), neighbors[points[i]])
    end

    return res
end


"""
_gen_graph(neighbors)

Given an array where res[i] <- the neighbors of node i, return the
corresponding lightgraph.
"""
function _gen_graph(neighbors)    
    res = LightGraphs.SimpleGraph(length(neighbors))

    for i in 1:length(neighbors)
        for j in 1:length(neighbors[i])
            LightGraphs.add_edge!(res, i => neighbors[i][j])
        end
    end

    return res
end

"""
    plot_field(field)

Plots the obstacles to the existing active plot.
"""
function plot_field(field)
    for i = 1:length(field)
        Plots.plot!(field[i], xlims = (-0.05,1.05), ylim = (-0.05, 1.05))
    end
end

"""
    plot_lines(field)

Plots the lines that make up the obstacles' halfspaces to the existing active
plot.
"""
function plot_lines(field)
    halfspaces = @pipe map(obstacle -> Polyhedra.hrep(obstacle).halfspaces, field) |> Iterators.flatten(_) |> collect(_)
    unique!(halfspaces)

    for h in halfspaces
        if abs(h.a[2]) != 0//1
            f = x -> (h.β - h.a[1] * x) / h.a[2]

            # x = 0//1:1//100:1//1
            x = Rational.(LinRange(0,1,11))
            y = map(f, x)
        else
            x = [abs( h.β / h.a[1]); abs(h.β / h.a[1])]
            y = [0//1; 1//1]
        end

        Plots.plot!(x, y)
    end
end

"""
    plot_borders()

Plots the lines that make up the unit square's borders to the existing active
plot.
"""
function plot_borders()
    halfspaces = [
        Polyhedra.HalfSpace([1, 0], 0),
        Polyhedra.HalfSpace([1, 0], 1),
        Polyhedra.HalfSpace([0, 1], 0),
        Polyhedra.HalfSpace([0, 1], 1)
    ]

    for h in halfspaces
        f = x -> (h.β - h.a[1] * x) / h.a[2]

        x = 0:0.01:1
        y = map(f, x)

        Plots.plot!(x, y)
    end
end

"""
    plot_intersections(field)

Plots and numbers the points where the lines that make up the obstacles'
halfspaces to the existing active plot.
"""
function plot_intersections(field)
    intersections = find_intersections(field)[1]

    x = map(point -> point[1], intersections)
    y = map(point -> point[2], intersections)

    return Plots.scatter!(x,y, series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:length(x)]))
end

"""
    plot_new(n, name)

Generates a new obstacle course with n obstacles and plots it, saving it to
name.png
"""
function plot_new(n::Int, name::String, seed::Int64=11)
    Plots.plot()

    obs = gen_field(n,seed)

    plot_field(obs)
    plot_lines(obs)
    a = plot_intersections(obs)

    Plots.savefig(a, name)

    return construct_graph(obs)
end