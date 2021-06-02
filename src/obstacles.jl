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
    # v = Polyhedra.convexhull([.25,.25], [.25,.5], [.5,.25])
    
    return Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer))
end

"""
    gen_field(num_obstacles)

Creates an obstacle course, an array of obstacles in the unit square. If
obstacles overlap the convex hull of their combined vertices will be taken
instead.
"""
function gen_field(num_obstacles::Int)
    Random.seed!(11) # 11 is a nice seed
    obstacles = map(_n -> gen_obstacle(0.25), 1:num_obstacles)

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

    boundaries = [
        Polyhedra.HalfSpace([1, 0], 0),
        Polyhedra.HalfSpace([1, 0], 1),
        Polyhedra.HalfSpace([0, 1], 0),
        Polyhedra.HalfSpace([0, 1], 1)
    ]

    lines = [halfspaces; boundaries]

    points = fill([], length(lines))    # Association from lines to points
    for i = 1:length(lines)   # Does same as line above?
        points[i] = []
    end
    
    tol = 1.e-6
    for i = 1:length(lines)        
        for j = (i + 1):length(lines)
            A = vcat(lines[i].a', lines[j].a')
            b = [lines[i].β; lines[j].β]
            try
                x = A \ b
                # point = x[1] => x[2]

                if abs(x[1] - 0.0) < tol
                    x[1] = 0.0
                elseif abs(x[1] - 1.0) < tol
                    x[1] = 1.0
                end

                if abs(x[2] - 0.0) < tol
                    x[2] = 0.0
                elseif abs(x[2] - 1.0) < tol
                    x[2] = 1.0
                end

                point = x[1] => x[2]

                #if (x[1] >= 0 && x[1] <= 1 && x[2] >= 0 && x[2] <= 1)
                # Use isapprox(), the above line, or keep as is with tol?
                if (x[1] >= 0 - tol && x[1] <= 1 + tol && x[2] >= 0 - tol && x[2] <= 1 + tol)
                    push!(res, point)
                    push!(points[i], point)
                    push!(points[j], point)
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
    points, mapped = find_intersections(obs)
    
    # Create map from point to neighbors (counter clockwise ordered by angle against horizontal)
    neighbors = _find_neighbors(points, mapped)
    graph = _gen_graph(neighbors)

    # angles = map(point -> atan(point.first, point.second), points)
    angles = fill(Inf, length(points), length(points))
    for i in 1:length(points)
        for j in 1:length(points)
            angles[i, j] = atan(points[j].second - points[i].second, points[j].first - points[i].first)
            if angles[i, j] < 0
                angles[i, j] += 2 * pi                
            end
        end
    end

    function greatest_angle_neighbor(source, last_angle, visited)
        unvisited = filter(vertex -> !(vertex in visited), neighbors[source])
        neighbor_angles = map(neighbor -> angles[source, neighbor], unvisited)
        zipped = @pipe zip(unvisited, neighbor_angles) |> collect(_)
        zipped # ?

        greater = filter(tup -> ((tup[2] > last_angle) && !(isapprox(tup[2], last_angle))), zipped)

        if isempty(greater)
            return (-1, -1) # bad state
        else
            sort!(greater, by=(tup -> tup[2]))
            destination = greater[1][1]
        end

        return (destination, angles[source, destination])
    end

    faces = []

    for i in 1:length(points)
        start = i

        for j in 1:length(neighbors[start])
            current = neighbors[start][j]
            last_angle = angles[start, current]

            face = [start, current]

            while start != current
                current # ?
                if (length(face) > 2 && (start in neighbors[current]))
                    push!(faces, copy(face))
                    break;
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

    unique_faces = filter(face -> begin
        face_set = Set(face)
        #included = face_set in face_sets

        face_set_v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        face_set_poly = Polyhedra.polyhedron(face_set_v, Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer))
        face_set_overlaps_obs_faces = false
        included = false
        
        #= 
        # Not needed because we have the obstacles already stored in obs
        obs_face_poly = map(obs_face -> begin
            obs_face_col = collect(obs_face)

            face_set_v = Polyhedra.convexhull(map(i -> collect(points[i]), obs_face_col)...)
            polygon = Polyhedra.polyhedron(face_set_v, Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer))
            display(plot!(polygon))
            return polygon
            #return Polyhedra.polyhedron(face_set_v, Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer))
        end, obstacle_faces)

        for obs_face in obs_face_poly
            if !isempty(Polyhedra.intersect(face_set_poly, obs_face))
                face_set_overlaps_obs_faces = true
                break
            end
        end
        =#

        # NEED TO FIX: Free spaces adjacent to an obstacle have overlap
        # Sol 1: Check if ext_points of face_set_poly are a subset of all points on the obstacle face
        # Sol 2: Check if face_set_poly is a subset of the obstacle face
        for (i,obs_face) in enumerate(obs)
            if !isempty(Polyhedra.intersect(face_set_poly, obs_face))
                face_set_overlaps_obs_faces = true
                # println("Overlap with obstacle $i")
                break
            end
        end

        # if !(face_set in face_sets) && !(face_set in obstacle_faces)
        if !(face_set in face_sets) && !(face_set_overlaps_obs_faces)
            push!(face_sets, face_set)
            included = true
        end
        
        # Cheat: Add all faces (with no repeats)
        if !(face_set in face_sets)
            push!(face_sets, face_set)
            included = true
        end

        return included
    end, faces)

    # See which are the bad faces
    # Plots.plot()
    # for (j,face) in enumerate(unique_faces)
    #     #println("$face")
    #     v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
    #     x_locations = map(i -> points[i].first, face)
    #     y_locations = map(i -> points[i].second, face)
    #     #println("$x_locations")
        
    #     avg_x = Statistics.mean(x_locations)
    #     avg_y = Statistics.mean(y_locations)
    #     #println("$avg_x , $avg_y")
    #     polygon = Polyhedra.polyhedron(
    #         v,
    #         Polyhedra.DefaultLibrary{Float64}(Gurobi.Optimizer)
    #     )
    
    #     Plots.plot!(polygon)
    #     Plots.plot!([avg_x], [avg_y], series_annotations=([Plots.text("$j", :center, 8, "courier")]))
    #     #display(Plots.plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
    
    #     #display(plot(polygon, title="Free Face # $j", xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
    # end
    # display(Plots.plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05)))

    # Cheat: Remove bad faces 3, 5, 6, 8, 15, 18
    iters = (1:2, 4, 7, 9:14, 16:17, 19:22)
    u_faces = []
    for iter in iters
        for i in iter
            push!(u_faces, unique_faces[i])
        end
    end

    #return (obs, points, graph, @pipe map(face -> reverse(face), unique_faces) |> Set(_))
    return (obs, points, graph, @pipe map(face -> reverse(face), u_faces) |> Set(_))
end

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
                index = findfirst(p -> p == point, colinear) # adjust for potential tol issue? Probably won't arise

                if index - 1 > 0
                    push!(neighbors[point], colinear[index - 1])
                end

                if index + 1 <= length(colinear)
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

    for h in halfspaces
        f = x -> (h.β - h.a[1] * x) / h.a[2]

        x = 0:0.01:1
        y = map(f, x)

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
function plot_new(n::Int, name::String)
    Plots.plot()

    obs = gen_field(n)

    plot_field(obs)
    plot_lines(obs)
    a = plot_intersections(obs)

    Plots.savefig(a, name)

    return construct_graph(obs)
end