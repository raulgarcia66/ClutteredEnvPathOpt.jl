import LightGraphs
import Polyhedra
import GLPK
import Plots

using Pipe


function gen_obstacle()
    max_side_len = 0.25
    point = rand(2)

    v = Polyhedra.convexhull(point, map(x -> x > 1 ? 1 : x, point + max_side_len * rand(2)), map(x -> x > 1 ? 1 : x, point + max_side_len * rand(2)))
    return Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer))
end

function gen_field(num_obstacles::Int)
    obstacles = map(_n -> gen_obstacle(), 1:num_obstacles)

    return remove_overlaps(obstacles)
end

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

function has_overlaps(obstacles)
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

function find_intersections(obstacles)
    res = Set([])

    halfspaces = @pipe map(obstacle -> Polyhedra.hrep(obstacle).halfspaces, obstacles) |> Iterators.flatten(_) |> collect(_)

    boundries = [
        Polyhedra.HalfSpace([1, 0], 0),
        Polyhedra.HalfSpace([1, 0], 1),
        Polyhedra.HalfSpace([0, 1], 0),
        Polyhedra.HalfSpace([0, 1], 1)
    ]

    lines = [halfspaces; boundries]

    points = fill([], length(lines))    # Association from lines to points
    for i = 1:length(lines)
        points[i] = []
    end
    
    for i = 1:length(lines)        
        for j = (i + 1):length(lines)
            A = vcat(lines[i].a', lines[j].a')
            b = [lines[i].β; lines[j].β]
            try
                x = A \ b
                point = x[1] => x[2]

                if (x[1] >= 0 && x[1] <= 1 && x[2] >= 0 && x[2] <= 1)
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

function plot_field(field)
    for i = 1:length(field)
        Plots.plot!(field[i], ylim = (0, 1))
    end
end

function plot_lines(field)
    halfspaces = @pipe map(obstacle -> Polyhedra.hrep(obstacle).halfspaces, field) |> Iterators.flatten(_) |> collect(_)

    for h in halfspaces
        f = x -> (h.β - h.a[1] * x) / h.a[2]

        x = 0:0.01:1
        y = map(f, x)

        Plots.plot!(x, y)
    end
end

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

function plot_intersections(field)
    intersections = find_intersections(field)[1]

    x = map(point -> point[1], intersections)
    y = map(point -> point[2], intersections)

    # Plots.scatter!(x, y)
    Plots.scatter!(x,y, series_annotations=([Plots.text(string(x), :right, 6, "courier") for x in 1:length(x)]))
end

function plot_new(n)
    Plots.plot()

    obs = gen_field(n)

    plot_field(obs)
    plot_lines(obs)
    plot_intersections(obs)

    construct_graph(obs)
end

function construct_graph(obs)
    @show points, mapped = find_intersections(obs)
    
    # Create map from point to neighbors (counter clockwise ordered by angle against horizontal)
    @show neighbors = find_neighbors(points, mapped)
    @show graph = gen_graph(neighbors)

    # angles = map(point -> atan(point.first, point.second), points)
    angles = fill(Inf, length(points), length(points))
    for i in 1:length(points)
        for j in 1:length(points)
            angles[i, j] = atan(points[j].second - points[i].second, points[j].first - points[i].first)
            if angles[i, j] <= 0
                angles[i, j] += 2 * pi                
            end
        end
    end
    @show angles

    faces = []

    for i in 1:length(points)
        @show start = i

        for j in 1:length(neighbors[start])
            @show current = neighbors[start][j]
            @show last_angle = angles[start, current]

            bad = false
            face = [start, current]

            while start != current
                @show current
                if (length(face) > 2 && start in neighbors[current])
                    @show "POOP"
                    break;
                end

                @show viable = filter(neighbor -> !(neighbor in face) && angles[current, neighbor] > last_angle && !isapprox(angles[current, neighbor], last_angle), neighbors[current])
                # if (length(unvisited) == 1)
                #     @show last_angle = angles[current, unvisited[1]]
                #     @show current = unvisited[1]
                # else
                #     @show last_angle = @pipe map(vertex -> angles[current, vertex] > last_angle && !isapprox(angles[current, vertex], last_angle) ? angles[current, vertex] : Inf, unvisited) |>
                #                              min(_...)
                #     @show current = filter(vertex -> angles[current, vertex] == last_angle, unvisited)[1]

                # end

                temp = sort(viable, by=(vertex -> angles[current, vertex]))
                if isempty(temp)
                    bad = true
                    break
                end
                next = temp[end]
                last_angle = angles[current, next]
                current = next
                # neighbor_angles = map(neighbor -> angles[current, neighbor], unvisited)
                # viable_angles = filter(angle -> angle >= last_angle, neighbor_angles)
                # last_angle = min(viable_angles...)[1]
                # current = @pipe findall(angle -> angle == last_angle, angles) |> filter(vertex -> vertex in unvisited, _)[1]

                push!(face, current)
            end

            if !bad
                push!(faces, copy(face))
            end
            @show faces
        end
    end

    @show faces

    # Pick arbitrary node A
    # Pop its largest angle neighbor off mapped stack and go to it
    # Repeat until returned to A
    # Save face data
    # continue until A's stack is empty
    # Pick new arbitrary point A
    # Continue until no points remain
end

function find_neighbors(points, mapped)
    neighbors = Dict()

    for point in points
        neighbors[point] = []

        for colinear in mapped
            if in(point, colinear)
                index = findfirst(p -> p == point, colinear)

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

function gen_graph(neighbors)    
    res = LightGraphs.SimpleGraph(length(neighbors))

    for i in 1:length(neighbors)
        for j in 1:length(neighbors[i])
            LightGraphs.add_edge!(res, i => neighbors[i][j])
        end
    end

    return res
end

plot_new(1)