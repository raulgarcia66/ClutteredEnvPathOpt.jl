import LightGraphs
import Polyhedra
import GLPK
import Random
import Statistics
import Triangulate
# using Pipe
# import Plots

"""
    gen_obstacle(max_side_len)

Creates an obstacle from a given set of points in the unit square.
"""
function gen_obstacle(vertices...)
    # Points can be given in any order since they are ordered clockwise within the polyhedron generated

    if length(vertices) < 4
        error("Obstacles need to have at least 4 extreme points.")
    elseif typeof(vertices[1][1]) !== Rational{Int64}
        error("Obstacle vertices currently need to be of type Rational{Int64}.")
    end

    v = Polyhedra.convexhull(vertices...)
    poly = Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
    Polyhedra.removevredundancy!(poly)

    if Polyhedra.npoints(poly) < 4
        error("Obstacles need to have at least 4 extreme points. Currently at $(Polyhedra.npoints(poly))")
        # Autofix: Can add a point between two points that is slightly shifted and not in the convex hull
        # points = Polyhedra.points(poly)
        # new_point = [ ; ]
        # v = Polyhedra.convexhull(points, new_point)
        # poly = Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
        # Polyhedra.removevredundancy!(poly)
    end
    
    Polyhedra.hrep(poly)
    return poly
end

"""
    gen_obstacle_four(point_per_dim)

Creates an obstacle with four extreme points from a grid of points in the unit square.
"""
function gen_obstacle_four(points_per_dim::Int)
    points_per_dim = max(points_per_dim, 6)  # Need at least 6 interior points
    nums = LinRange(0//1, 1//1, points_per_dim+2)[2:end-1]  # Do not want obstacles touching the unit square boundaries
    
    # Select the bottom left point at random and select the other three points with respect to it
    x_index = rand(2:(length(nums)-4))
    y_index = rand(2:(length(nums)-4))
    bottom_left = [ nums[x_index] ; nums[y_index] ]
    bottom_right = [ nums[x_index + 3 + rand(-1:1)] ; nums[y_index + rand(-1:1)] ]
    top_left = [ nums[x_index + rand(-1:1)] ; nums[y_index + 3 + rand(-1:1)] ]
    top_right = [ nums[x_index + 3 + rand(-1:1)] ; nums[y_index + 3 + rand(-1:1)] ]

    v = Polyhedra.convexhull(top_left, top_right, bottom_left, bottom_right)
    poly = Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
    Polyhedra.removevredundancy!(poly)
    # Loop below assures the polyhedron has 4 extreme points (sometimes the 4th is a convex comb.) for pairwise IB purposes
    while Polyhedra.npoints(poly) < 4
        top_right = [ nums[x_index + 3 + rand(-1:1)] ; nums[y_index + 3 + rand(-1:1)] ]
        v = Polyhedra.convexhull(top_left, top_right, bottom_left, bottom_right)
        poly = Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
        Polyhedra.removevredundancy!(poly)
    end

    Polyhedra.hrep(poly)
    return poly
end

"""
    gen_field_random(num_obstacles; custom, seed)

Creates an obstacle course, an array of obstacles in the unit square. If
obstacles overlap, the convex hull of their combined extreme points will be taken
instead.
"""
function gen_field_random(num_obstacles::Int; points_per_dim::Int = 20, seed::Int=1) where {T}
    Random.seed!(seed)
    obstacles = map(_n -> gen_obstacle_four(points_per_dim), 1:num_obstacles)

    return remove_overlaps(obstacles)
end

"""
    gen_field(obstacles)

Creates an obstacles course from a collection of obstacles by removing overlaps.
"""
function gen_field(obstacles)
    # TODO: May just remove this function
    remove_overlaps(obstacles)
end

"""
    remove_overlaps(obstacles)

Given an obstacle course, check for overlapping obstacles and replace offending
obstacles with the convex hull of their combined extreme points.
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
                    Polyhedra.removevredundancy!(accumulator)
                    Polyhedra.removehredundancy!(accumulator)
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

Analyzes the halfspaces that make up each obstacle and unit square border, and finds
where they intersect. Returns an array of all intersections, as well as an array
mapping each halfspace to the points of intersection that lie on it, sorted. Points
in the interior of obstacles are stored at the end of the array.
"""
function find_intersections(obstacles)
    res = Pair{Rational{Int64},Rational{Int64}}[]
    
    halfspaces = @pipe map(obstacle -> Polyhedra.hrep(obstacle).halfspaces, obstacles) |> Iterators.flatten(_) |> collect(_)
    unique!(halfspaces)

    boundaries = [
        Polyhedra.HalfSpace([-1//1, 0//1], 0//1),
        Polyhedra.HalfSpace([1//1, 0//1], 1//1),
        Polyhedra.HalfSpace([0//1, -1//1], 0//1),
        Polyhedra.HalfSpace([0//1, 1//1], 1//1)
    ]
    filter!(b -> !(b in halfspaces), boundaries)

    lines = [halfspaces; boundaries]

    # Save indices of horizontal and vertical halfspaces to keep track of which points need
    # the same y- or x- coordinates (necessary to sort correctly)
    vert_ind = []
    hor_ind = []
    for i = 1:length(lines)
        if abs(lines[i].a'[1]) == 1//1 && abs(lines[i].a'[2]) == 0//1 # set tolerance?
            push!(vert_ind, i)
        end
        if abs(lines[i].a'[1]) == 0//1 && abs(lines[i].a'[2]) == 1//1
            push!(hor_ind, i)
        end
    end

    points = fill(Pair{Rational{Int64},Rational{Int64}}[], length(lines))    # For storing the points on each halfspace
    for i = 1:length(lines)
        points[i] = Pair{Rational{Int64},Rational{Int64}}[]
    end

    # Function for determining if a point should be considered a duplicate of another. If so, the index
    # of the existing one is returned
    function point_duplicates(point, res)
        tol = Rational(1.e-13)
        for index in 1:length(res)
            if abs(point.first - res[index].first) < tol && abs(point.second - res[index].second) < tol 
                return true, index
            end
        end
        return false, -1
    end
    
    tol = Rational(1.e-13)
    # keys are the indices in "lines"
    x_vert = Dict([])
    y_hor = Dict([]) 
    for i = 1:length(lines)     
        for j = (i + 1):length(lines)
            A = vcat(lines[i].a', lines[j].a')
            b = [lines[i].β; lines[j].β]
            try
                x = A \ b

                # Assures boundary points do not have floating point errors
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

                if (x[1] >= 0//1 && x[1] <= 1//1 && x[2] >= 0//1 && x[2] <= 1//1)
                    # Assures points on vertical or horizontal halfspaces have the same x- or y- coordinates, resp.
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

                    duplicate, index = point_duplicates(point, res)
                    if !(duplicate)
                        # If point is not a duplicate, push it into res and the halfspaces it is in
                        push!(res, point)
                        push!(points[i], point)
                        push!(points[j], point)
                    else
                        # If point is a duplicate, push it into the halfspaces it's in if not added yet
                        if !(res[index] in points[i])
                            push!(points[i], res[index])
                        end
                        if !(res[index] in points[j])
                            push!(points[j], res[index])
                        end
                    end
                end
            catch e
                # No intersections found
            end
        end
    end

    # Find vertices in the interior of an obstacle
    points_vec = res
    inside_indices = filter(i -> begin
        inside = false
        for ob in obstacles
            if Polyhedra.ininterior([points_vec[i].first ; points_vec[i].second], ob)
                inside = true
                # println("Inside Point: $(points_vec[i].first), $(points_vec[i].second)")
                break
            end
        end
        return inside
    end, 1:length(points_vec))

    # Move inside vertices to the end of the points vector, so that they may be removed from the points
    # vector and the lightgraph graph without altering the vertex labels and the edges remaining. 
    # Removal does not occur here.
    for i = 1:(length(inside_indices))
        temp = points_vec[inside_indices[i]]
        points_vec[inside_indices[i]] = points_vec[end - (i-1)]
        points_vec[end - (i-1)] = temp
    end

    sorted_points = map(list -> sort(list), points)

    # Third return value is for knowing how many inside points there are
    return (points_vec, sorted_points, length(inside_indices))
end

"""
    find_points(obs)
Compute the points corresponding to obstacle extreme points and the four corners
of the unit square. Returns the points, and possibly the obstacle face vertices and edges.
Function assumes obstacles have rational data.
"""
function find_points(obs; with_faces::Bool=false)
    S = Vector{Int64}[]   # Store segments
    obstacle_faces = Vector{Int64}[]   # Store obstacle face vertices (should be ordered)
    points = Pair{Rational{Int64},Rational{Int64}}[]   # Store points as Vector of Pairs
    vertex_num = 0
    for ob in obs
        vertices_in_ob = Int64[]
        iter = Polyhedra.points(ob)
        # Store x and y coordinates in points
        for point in iter
            push!(points, point[1] => point[2])
            vertex_num += 1
            push!(vertices_in_ob, vertex_num)
        end
        push!(obstacle_faces, vertices_in_ob)

        # Push edges of obstacle boundary
        # Vertices are returned clockwise when obtained via Polyhedra.points(ob)
        if with_faces
            for i = 1:(length(vertices_in_ob)-1)
                push!(S, [vertices_in_ob[i] ; vertices_in_ob[i+1]])
            end
            push!(S, [vertices_in_ob[end] ; vertices_in_ob[1]])
        end
    end

    points = [points; 0//1 => 0//1 ; 0//1 => 1//1; 1//1 => 1//1; 1//1 => 0//1]
    # TODO: Function relies on the points being unique
    # unique!(points)  # careful with doing this; are there any dependencies on point order?

    if with_faces
        # Add unit cell boundary segments
        num_points = length(points)
        for i = (num_points-3):(num_points-1)
            push!(S, [i ; i+1])
        end
        push!(S, [num_points ; num_points-3])

        return points, obstacle_faces, S
    else 
        return points
    end

end

"""
    construct_graph(obs)

Given a list of obstacles with no overlaps, finds the intersections of the lines
that make up each obstacle's halfspaces and constructs a graph with the
intersections as nodes and the halfspaces as edges. Returns a tuple containing
the obstacles, a list of the location of the nodes' coordinates, the graph, a set of
every face of obstacle space, and a set of every face of free space. (Each face is 
a clockwise-ordered vector of vertices).
"""
function construct_graph(obs)
    points, mapped, inside_quant = ClutteredEnvPathOpt.find_intersections(obs)

    # Create map from point to its neighbors
    neighbors = ClutteredEnvPathOpt._find_neighbors(points, mapped)
    # Create lightgraph
    graph = ClutteredEnvPathOpt._gen_graph(neighbors)

    # TODO: Angles become floats once campared with the irrational pi
    # angles[i,j] is the angle point j makes with point i
    angles = fill(Rational(Inf), length(points), length(points))
    for i in 1:length(points)
        for j in 1:length(points)
            angles[i, j] = Rational(atan(points[j].second - points[i].second, points[j].first - points[i].first))
            if angles[i, j] < 0//1
                angles[i, j] += 2 * pi # Rational(2*3.141592653589793)
            end
        end
    end

    function greatest_angle_neighbor(source, last_angle, visited)
        #tol = Rational(1.e-15)
        unvisited = filter(vertex -> !(vertex in visited), neighbors[source])
        neighbor_angles = map(neighbor -> angles[source, neighbor], unvisited)
        # For sorting angles
        for i in 1:length(neighbor_angles)
            if neighbor_angles[i] >= Rational(3.141592653589793)  # Use pi?
                neighbor_angles[i] -= Rational(2*3.141592653589793)
            end
        end
        zipped = @pipe zip(unvisited, neighbor_angles) |> collect(_)
        # zipped is an array of tuples of the form: (unvisited node, angle it makes with source)

        # Two cases to assure we are constructing the most compact face (ie, no face with subfaces)
        if last_angle >= 0//1 && last_angle < Rational(3.141592653589793)
            greater = filter(tup -> (tup[2] > last_angle) && (tup[2] < last_angle + Rational(3.141592653589793)), zipped)
        else
            greater = filter(tup -> (tup[2] > last_angle - Rational(2*3.141592653589793)) && (tup[2] < last_angle - Rational(3.141592653589793)), zipped)
        end

        if isempty(greater)
            return (-1, -1) # bad state
        else
            sort!(greater, by=(tup -> tup[2]))
            destination = greater[end][1]
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

    # Find irredundant faces
    face_sets = []
    unique_faces = filter(face -> begin
        face_set = Set(face)
        included = false

        if !(face_set in face_sets)
            push!(face_sets, face_set)
            included = true
        end

        return included
    end, faces)

    # Find free faces (by determining which do not overlap obstacles)
    free_faces = filter(face -> begin
        face_set = Set(face)   # face is already a set, luckily taking a Set of a Set doesn't do anything
        # TODO: It's not even used! Comment it and test, then remove

        face_set_v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        face_set_poly = Polyhedra.polyhedron(face_set_v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
        face_set_overlaps_obs_faces = false
        included = false

        # If the intersection of face and an obstacle is the face itself, than that face is a subface of an obstacle
        tol =  Rational(1.e-10)
        for n in 1:length(obs)  # need to loop in this manner for 'continue' to work properly
            indicator = 0
            intersect_poly = Polyhedra.intersect(face_set_poly, obs[n])
            Polyhedra.vrep(intersect_poly)  # compute V-representation

            if isempty(Polyhedra.points(intersect_poly)) == true
                continue  # intersection is empty
            end

            ext_p_intersect_poly = sort(collect(intersect_poly.vrep.points.points))
            ext_p_face_set_poly = sort(collect(face_set_poly.vrep.points.points))
            if length(ext_p_intersect_poly) != length(ext_p_face_set_poly)
                continue  # intersection is not the face
            end

            for i in 1:length(ext_p_face_set_poly)
                for j in 1:length(ext_p_face_set_poly[i])  # should be 1:2 (x and y coordinates)
                    if abs(ext_p_intersect_poly[i][j] - ext_p_face_set_poly[i][j]) > tol
                        indicator += 1  # the extreme points do not coincide at the x or y coordinate
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

        if !(face_set_overlaps_obs_faces)
            included = true
        end

        return included
    end, unique_faces)


    # Swap obstacle (sub)faces with whole faces, then remove inside nodes and inside edges
    # Chad Subset Approach

    # Remove inside vertices
    if inside_quant > 0
        length_p = length(points)
        for i = 1:inside_quant
            LightGraphs.rem_vertex!(graph, length_p - (i-1))
            pop!(points)
        end
    end

    # Determine the vertices that make up each obstacle's extreme points
    obs_euclid_points = map(ob -> Polyhedra.points(ob), obs)
    obstacle_faces = map(ob_points -> begin
        tol = Rational(1e-13)
        face_points = []
        for ob_point in ob_points
            vertex_index = findfirst(p -> abs(ob_point[1] - p.first) < tol && abs(ob_point[2] - p.second) < tol, points)
            if vertex_index !== nothing
                push!(face_points, vertex_index)
            end
        end
        return face_points
    end, obs_euclid_points)

    # Find all vertices associated with an obstacle face and remove inside edges

    boundary_edges = fill([], length(obs))  # store boundary edges of each obstacle
    for i = 1:length(obs)
        boundary_edges[i] = []
    end
    
    edges_to_remove = []  # edges crossing through an obstacle are to be removed
    for edge in LightGraphs.edges(graph)
        on_boundary = true
        for i = 1:length(obs)
            if Polyhedra.in([points[edge.src].first; points[edge.src].second], obs[i]) && Polyhedra.in([points[edge.dst].first; points[edge.dst].second], obs[i])
                # Edge is contained in obstacle i
                # println("Found Edge: $(edge.src) => $(edge.dst) with i = $i")
                push!(obstacle_faces[i], edge.src)
                push!(obstacle_faces[i], edge.dst)
                unique!(obstacle_faces[i])
                if Polyhedra.ininterior([(points[edge.src].first + points[edge.dst].first)//2 ; (points[edge.src].second + points[edge.dst].second)//2], obs[i])
                    # Edge crosses through obstacle i
                    #println("Removing edge $(edge.src) => $(edge.dst)")
                    push!(edges_to_remove, (edge.src, edge.dst))
                    on_boundary = false
                end
                if on_boundary   # this can be change to an if-else with the above
                    # println("Adding edge $(edge.src) => $(edge.dst) with i = $i")
                    push!(boundary_edges[i], edge)
                end
                break
            end
        end
    end

    for (a,b) in edges_to_remove
        # println("Removing edge $a => $b")
        LightGraphs.rem_edge!(graph, a, b)
    end

    # for i in 1:length(boundary_edges)
    #     println("Obstacle $i")
    #     for edge in boundary_edges[i]
    #         println("$edge")
    #     end
    # end

    # Now need to order points on obstacle faces counter-clockwise (they're then returned clockwise)

    # Start with left-most point (if multiple, use bottom-most)
    starting_points = []
    for ob_face in obstacle_faces
        x_coor = map(vertex -> points[vertex].first, ob_face)
        min_x = min(x_coor...)
        candidates = filter(vertex -> points[vertex].first == min_x, ob_face)
        if length(candidates) > 1
            y_coor = map(vertex -> points[vertex].second, candidates)
            min_y = min(y_coor...)
            candidates = filter(vertex -> points[vertex].second == min_y, candidates)
        end
        push!(starting_points, candidates[1])
    end

    # println("Starting Points:")
    # for sp in starting_points
    #     println("$sp")
    # end

    obstacle_faces_ordered = []
    for i = 1:length(starting_points)
        sp = starting_points[i]
        sp_adj = []  #  adjacent nodes (should be exactly 2)
        for edge in boundary_edges[i]
            if sp == edge.src
                push!(sp_adj, edge.dst)
            elseif sp == edge.dst
                push!(sp_adj, edge.src)
            end
        end
        sp_n1 = sp_adj[1]
        sp_n2 = sp_adj[2]

        angle_n1 = angles[sp, sp_n1]
        if angle_n1 >= Rational(3.141592653589793)  # for ordering
            angle_n1 -= Rational(2*3.141592653589793)
        end
        angle_n2 = angles[sp, sp_n2]
        if angle_n2 >= Rational(3.141592653589793)
            angle_n2 -= Rational(2*3.141592653589793)
        end

        # Determine counter-clockwise order of the three points
        n1 = sp
        if angle_n1 > angle_n2
            face = [sp_n1; sp; sp_n2] 
            n0 = sp_n1
            n2 = sp_n2
        else
            face = [sp_n2; sp; sp_n1]
            n0 = sp_n2
            n2 = sp_n1
        end

        # From current node n2, determine the next node in the face by checking which other node n2
        # is adjacent to (other than n1). Repeat until the next node is n0.
        while n2 != n0
            next_found = false
            for edge in boundary_edges[i]
                if edge.src == n2 && edge.dst != n1
                    n1 = copy(n2) 
                    n2 = edge.dst
                    next_found = true
                    break
                elseif edge.dst == n2 && edge.src != n1
                    n1 = copy(n2) 
                    n2 = edge.src
                    next_found = true
                    break
                end
            end
            if !next_found
                println("Infinite while loop. Face may not be complete.")
                break
            end
            if n2 != n0
                push!(face, n2)
            end
        end

        push!(obstacle_faces_ordered, face)
    end

    return (obs, points, graph, Set(map(face -> reverse(face), obstacle_faces_ordered)), Set(map(f -> reverse(f), free_faces)))
end

"""
    _find_neighbors(points, mapped)

Given a list of points and an association of halfspace lines to points that lie
on the line, return an array where res[i] <- the neighbors of node i.
"""
function _find_neighbors(points::Vector{Pair{Rational{Int64},Rational{Int64}}}, mapped::Vector{Vector{Pair{Rational{Int64},Rational{Int64}}}})::Vector{Vector{Int64}}
    # This function should only be used in construct_graph(), before "inside" points are removed
    neighbors = Dict()

    for point in points
        neighbors[point] = Pair{Rational{Int64},Rational{Int64}}[]

        for colinear in mapped
            if in(point, colinear)
                index = findfirst(p -> p == point, colinear)
                # Do not push the same point more than once
                if index - 1 > 0 && !(colinear[index - 1] in neighbors[point])
                    push!(neighbors[point], colinear[index - 1])
                end

                if index + 1 <= length(colinear) && !(colinear[index + 1] in neighbors[point])
                    push!(neighbors[point], colinear[index + 1])
                end
            end
        end
    end

    # Res will contain the neighbors by their node label
    res = fill(Int64[], length(points))
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
    construct_graph_delaunay(obs)

Generate the constrained Delaunay triangulation partition of the free space.
Returns a tuple containing the obstacles, a list of the location of the nodes' coordinates,
the graph, a set of every face of obstacle space, and a set of every face of free space.
(Each face is a clockwise-ordered vector of vertices).
"""
function construct_graph_delaunay(obs; merge_faces=true)
    # Code below was moved to find_points function
    # S = []   # Store segments
    # obstacle_faces = []   # Store obstacle face vertices (should be ordered)
    # points = []   # Store points as Vector of Pairs
    # vertex_num = 0
    # for ob in obs
    #     vertices_in_ob = []
    #     iter = Polyhedra.points(ob)
    #     # Store x and y coordinates in points
    #     for point in iter
    #         push!(points, point[1] => point[2])
    #         vertex_num += 1
    #         push!(vertices_in_ob, vertex_num)
    #     end
    #     # Push edges of obstacle boundary (assumes they're in order) 
    #     for i = 1:(length(vertices_in_ob)-1)
    #         push!(S, [vertices_in_ob[i] ; vertices_in_ob[i+1]])
    #     end
    #     push!(S, [vertices_in_ob[end] ; vertices_in_ob[1]])

    #     push!(obstacle_faces, vertices_in_ob)
    # end

    # points = [points; 0//1 => 0//1 ; 0//1 => 1//1; 1//1 => 1//1; 1//1 => 0//1]
    # unique!(points)   # Code above relies on the points being unique

    # # Add boundary segments
    # num_points = length(points)
    # for i = (num_points-3):(num_points-1)
    #     push!(S, [i ; i+1])
    # end
    # push!(S, [num_points ; num_points-3])

    points, obstacle_faces, S = ClutteredEnvPathOpt.find_points(obs; with_faces=true)

    # Initialize triangulate structure
    triin = Triangulate.TriangulateIO()

    # Store points in matrix
    C = zeros(2,length(points))
    for j = 1:length(points)
        C[:,j] = [points[j].first ; points[j].second]
    end
    triin.pointlist = Matrix{Cdouble}(C)

    # Store segments in matrix
    segments = zeros(Int,2,length(S))
    for j=1:length(S)
        segments[:,j] = [S[j][1]; S[j][2]]
    end
    triin.segmentlist = Matrix{Cint}(segments)

    # Add hole list
    hole_center = zeros(2,length(obstacle_faces))
    for (j,face) in enumerate(obstacle_faces)
        x_locations = map(i -> points[i].first, face)
        y_locations = map(i -> points[i].second, face)
        avg_x = Statistics.mean(x_locations)
        avg_y = Statistics.mean(y_locations)
        hole_center[:,j] = [avg_x; avg_y]
    end
    triin.holelist = Matrix{Cdouble}(hole_center)

    # Compute CDT
    (triout, _) = Triangulate.triangulate("pQ", triin)

    triangles = Matrix{Int64}(triout.trianglelist)   # matrix with 3 rows
    # edges = Matrix{Int64}(triout.edgelist)   # matrix with 2 rows

    # Create graph
    graph = LightGraphs.SimpleGraph(length(points))
    # Add edges; Use output matrix of edges instead?
    for j = 1:size(triangles,2)
        LightGraphs.add_edge!(graph, triangles[1,j] => triangles[2,j])
        LightGraphs.add_edge!(graph, triangles[1,j] => triangles[3,j])
        LightGraphs.add_edge!(graph, triangles[2,j] => triangles[3,j])
    end

    # Store free faces as vectors (utilizing holes ⟹ all triangles are free faces)
    free_faces = Vector{Int64}[]
    for j = 1:size(triangles,2)
        push!(free_faces, reverse(triangles[:,j]))
    end

    # Face merging
    if merge_faces
        free_faces, graph = face_merger(free_faces, graph, points)
    end

    # return (obs, points, graph, Set(map(face -> face, obstacle_faces)), Set(map(f -> reverse(f), free_faces)))
    return (obs, points, graph, Set(map(face -> face, obstacle_faces)), Set(map(f -> f, free_faces)))
end

"""
    face_merger(faces, graph, points)
Given a set of vectors corresponding to faces, this function merges faces together when possible
to decrease the number of faces and edges in the skeleton graph. 
"""
function face_merger(faces::Vector{Vector{T}}, graph::W, points::Vector{Pair{Rational{Int64},Rational{Int64}}}) where {T, W <: LightGraphs.AbstractGraph}
    edges_to_remove = Tuple{T,T}[]
    
    merge_found = true
    while merge_found
        sort!(faces, by=length)
        merge_found = false
        for i = 1:length(faces)
            for j = i+1:length(faces)
                # merge_found = false
                face_i = Set(faces[i])
                face_j = Set(faces[j])
                if length(intersect(face_i, face_j)) == 2
                    v_i = Polyhedra.convexhull(map(i -> collect(points[i]), faces[i])...)
                    v_j = Polyhedra.convexhull(map(i -> collect(points[i]), faces[j])...)
                    new_face_vec = union(face_i, face_j)   # not really a vector
                    v_new = Polyhedra.convexhull(map(i -> collect(points[i]), collect(new_face_vec))...)

                    polygon_i = Polyhedra.polyhedron(
                        v_i,
                        Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer)
                    )
                    polygon_j = Polyhedra.polyhedron(
                        v_j,
                        Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer)
                    )
                    polygon_new = Polyhedra.polyhedron(
                        v_new,
                        Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer)
                    )
                    if (Polyhedra.volume(polygon_i) + Polyhedra.volume(polygon_j)) ≈ Polyhedra.volume(polygon_new)
                        Polyhedra.removevredundancy!(polygon_new)
                        if Polyhedra.npoints(polygon_i) + Polyhedra.npoints(polygon_j) - 2 == Polyhedra.npoints(polygon_new)
                            merge_found = true
                            # Find new face
                            shared_edge_vec = intersect(faces[i], faces[j])
                            push!(edges_to_remove, (shared_edge_vec[1], shared_edge_vec[2]))
                            if shared_edge_vec == vcat(faces[i][1], faces[i][end]) && length(faces[i]) > 2
                                faces[i] = circshift(faces[i],1)
                                reverse!(shared_edge_vec)   # makes shared edge same order as faces[i]
                            end
                            first_index = findfirst(vertex -> vertex == shared_edge_vec[1],faces[i])
                            second_index = first_index + 1
                            # Shift faces[j] so that the second vertex is at end (first vertex is then at beginning)
                            index_of_second_vertex_in_edge = findfirst(v -> v == shared_edge_vec[2], faces[j])
                            shift_size = length(faces[j]) - index_of_second_vertex_in_edge
                            faces[j] = circshift(faces[j], shift_size)

                            new_face_vec_final = vcat(faces[i][1:first_index], faces[j][2:end-1], faces[i][second_index:end])
                            push!(faces, new_face_vec_final)

                            deleteat!(faces,j)
                            deleteat!(faces,i)

                            break
                        end
                    end
                end
            end
            if merge_found == true
                break
            end
        end
    end

    # Remove edges from graph
    for (a,b) in edges_to_remove
        LightGraphs.rem_edge!(graph, a, b)
    end

    return faces, graph
end