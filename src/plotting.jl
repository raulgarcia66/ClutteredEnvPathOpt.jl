import Plots
import Statistics
import LightGraphs
import Polyhedra
import GLPK

"""
    plot_faces(faces, points; plot_name, col, new_plot, individually)

Given a set of vectors corresponding to faces of a planar graph, plots the faces in the unit square.
"""
function plot_faces(faces::Set{Vector{T}}, points::Vector{Pair{Rational{Int64},Rational{Int64}}}; plot_name::String="Faces", col::String="green3", new_plot::Bool=true, individually::Bool=false) where {T}
    if new_plot && !individually
        Plots.plot()
    end

    for (j,face) in enumerate(faces)
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        x_locations = map(i -> points[i].first, face)
        y_locations = map(i -> points[i].second, face)

        # avg_x = Statistics.mean(x_locations)
        # avg_y = Statistics.mean(y_locations)

        polygon = Polyhedra.polyhedron(
            v,
            Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer)
        )

        if individually
            # Plot each face on a separate plot with node labels
            Plots.plot(polygon, color=col, alpha=0.9)
            Plots.plot!(x_locations, y_locations, series_annotations=([Plots.text(string(x), :center, 8, "courier") for x in face]))
            display(Plots.plot!(xlims=(-0.05,1.05), ylims=(-0.05,1.05), title="Face $j: $face"))
        else
            # Accumulate the faces on the same plot
            Plots.plot!(polygon, color=col, alpha=0.9)
            # Plots.plot!([avg_x], [avg_y], series_annotations=([Plots.text("$j", :center, 8, "courier")]))
            # display(Plots.plot!(polygon, title="Free Face # $j", xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
        end
    end

    if !individually
        display(Plots.plot!(title=plot_name, xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
    end
end

"""
    plot_edges(lg, points; plot_name, col, new_plot)

Given a LabeledGraph, plots its nodes and edges in the unit square.
"""
function plot_edges(lg::LabeledGraph{T}, points::Vector{Pair{Rational{Int64},Rational{Int64}}}; plot_name::String="Edges", col::String="colorful", new_plot::Bool=true, vertices::Dict{T,T}=Dict{T,T}()) where {T}
    if new_plot
        Plots.plot()
    end

    # Need to map nodes to the point they refer to in "points"
    rev = ClutteredEnvPathOpt._reverse_labels(lg.labels)
    for edge in LightGraphs.edges(lg.graph)
        if col == "colorful"
            Plots.plot!([points[rev[edge.src]].first, points[rev[edge.dst]].first], [points[rev[edge.src]].second, points[rev[edge.dst]].second],linewidth=2)
            # display(Plots.plot!(title="Edge ($(rev[edge.src]), $(rev[edge.dst]))"))
        else
            Plots.plot!([points[rev[edge.src]].first, points[rev[edge.dst]].first], [points[rev[edge.src]].second, points[rev[edge.dst]].second], color=col,linewidth=2)
            # display(Plots.plot!(title="Edge ($(rev[edge.src]), $(rev[edge.dst]))"))
        end
    end

    # ClutteredEnvPathOpt.plot_intersections(obstacles, vertices=vertices)
    ClutteredEnvPathOpt.plot_points(points, vertices=vertices)

    display(Plots.plot!(title=plot_name, xlims=(-0.05,1.05), ylims=(-0.05,1.05), legend=false))
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
            f = x -> (h.β - h.a[1] * x) // h.a[2]

            x = LinRange(0//1, 1//1, 11)
            y = map(f, x)
        else  # vertical lines
            x = fill(abs(h.β // h.a[1]), 11)
            y = LinRange(0//1, 1//1, 11)
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
        Polyhedra.HalfSpace([-1//1, 0//1], 0//1),
        Polyhedra.HalfSpace([1//1, 0//1], 1//1),
        Polyhedra.HalfSpace([0//1, -1//1], 0//1),
        Polyhedra.HalfSpace([0//1, 1//1], 1//1)
    ]

    for h in halfspaces
        if abs(h.a[2]) != 0//1
            f = x -> (h.β - h.a[1] * x) // h.a[2]

            x = LinRange(0//1, 1//1, 11)
            y = map(f, x)
        else  # vertical lines
            x = fill(abs(h.β // h.a[1]), 11)
            y = LinRange(0//1, 1//1, 11)
        end

        Plots.plot!(x, y)
    end
end

"""
    plot_intersections(field)

Plots and labels the points where the lines that make up the obstacles'
halfspaces intersect to the existing active plot.
"""
function plot_intersections(field; vertices::Dict{T,T}=Dict{T,T}()) where {T}
    intersections, _, inside_quant = find_intersections(field)

    # Remove points inside obstacles
    for _ = 1:inside_quant
        pop!(intersections)
    end

    if !isempty(vertices)
        points_filtered = []
        for v in keys(vertices)
            push!(points_filtered, intersections[v])
        end
        x = map(point -> point.first, points_filtered)
        y = map(point -> point.second, points_filtered)
        Plots.scatter!(x,y, color="red", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in keys(vertices)]))
    else
        x = map(point -> point.first, intersections)
        y = map(point -> point.second, intersections)
        Plots.scatter!(x,y, color="red", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:length(points)]))
    end

    # TODO: Delete this once certain it works
    # if !isempty(vertices)
    #     intersections = @pipe filter(v -> v in keys(vertices), 1:length(intersections)) |> intersections[_]
    # end

    # x = map(point -> point[1], intersections)
    # y = map(point -> point[2], intersections)

    # # Plots.scatter!(x,y, color="red3", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:length(x)]))
    # Plots.scatter!(x,y, color="red", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:length(x)]))
end

"""
    plot_points(points; vertices)

Plots and labels points given.
"""
function plot_points(points::Vector{Pair{Rational{Int64},Rational{Int64}}}; vertices::Dict{Int,Int}=Dict{Int,Int}())
    
    Plots.plot!(legend=false)

    if !isempty(vertices)
        points_filtered = []
        for v in keys(vertices)
            push!(points_filtered, points[v])
        end
        x = map(point -> point.first, points_filtered)
        y = map(point -> point.second, points_filtered)
        Plots.scatter!(x,y, color="red", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in keys(vertices)]))
    else
        x = map(point -> point.first, points)
        y = map(point -> point.second, points)
        Plots.scatter!(x,y, color="red", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:length(points)]))
    end

end

"""
    plot_new(n, name; custom, seed, save_image)

Generates a new obstacle course with n obstacles. If custom = true, obstacles will be
created from an array containing the points for each obstacle. If save_image = true,
the course will be plotted and saved to an image.
"""
function plot_new(n::Int, name::String; custom::Bool=false, seed::Int=1, save_image::Bool=false, partition::String="CDT", merge_faces=true)
    if custom
        # obs = gen_field(num_obstacles::Int; custom::Bool=false, points_it::Set{Vector{T}}=Set{Vector{Rational}}(), seed::Int=1) where {T}
    else
        obs = gen_field(n, seed = seed)
    end

    if save_image
        Plots.plot()
        plot_field(obs)
        points = ClutteredEnvPathOpt.find_points(obs)
        ClutteredEnvPathOpt.plot_points(points)
        #plot_lines(obs)
        #plot_borders()
        #plot_intersections(obs)
        Plots.png(name)
        display(Plots.plot!(title="Field"))
    end

    if partition == "CDT"
        return construct_graph_delaunay(obs; merge_faces=merge_faces)
    else
        return construct_graph(obs)
    end
end
