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
function plot_edges(lg::LabeledGraph{T}, points::Vector{Pair{Rational{Int64},Rational{Int64}}}; plot_name::String="Edges", col::String="colorful", new_plot::Bool=true, vertices::Dict{T,T}=Dict{T,T}(), with_labels=true) where {T}
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
    ClutteredEnvPathOpt.plot_points(points, vertices=vertices, with_labels=with_labels)

    display(Plots.plot!(title=plot_name, xlims=(-0.05,1.05), ylims=(-0.05,1.05), legend=false))
end

"""
    plot_biclique_cover(lg, points, cover; with_all)

Given a LabeledGraph and a biclique cover, plot each biclique.
"""
function plot_biclique_cover(lg::LabeledGraph{T}, points::Vector{Pair{Rational{Int64},Rational{Int64}}}, cover::Set{Pair{Set{T}, Set{T}}}; with_all::Bool=false, name::String="Biclique") where {T}
    e_bar = LightGraphs.edges(LightGraphs.complement(lg.graph))

    # Vector of sets of pairs (edges)
    cover_vec = collect(cover)
    BC_edges = map(pair -> ClutteredEnvPathOpt._cartesian_product(pair.first, pair.second), cover_vec)
    # all_BC_edges = reduce(union!, BC_edges, init=Set{Pair{T,T}}())
    # temp_graph = LightGraphs.SimpleGraph(LightGraphs.nv(lg.graph))   # will contain BC edges
    rev = ClutteredEnvPathOpt._reverse_labels(lg.labels)

    # for edge in all_BC_edges
    #     # edge.first/second will be our nodes, so we map them to the nodes they are in the graph
    #     LightGraphs.add_edge!(temp_graph, lg.labels[edge.first], lg.labels[edge.second])
    # end

    colors = [:firebrick1, :dodgerblue, :limegreen]

    # Plot biclique
    for (j,biclique_edges) in enumerate(BC_edges)
        plot(legend=false)
        if with_all
            temp_graph = LightGraphs.SimpleGraph(LightGraphs.nv(lg.graph))
            for edge in biclique_edges
                # edge.first/second will be our nodes, so we map them to the nodes they are in the graph
                LightGraphs.add_edge!(temp_graph, lg.labels[edge.first], lg.labels[edge.second])
            end

            for edge in e_bar
                if !(edge in LightGraphs.edges(temp_graph))
                    Plots.plot!([points[rev[edge.src]].first, points[rev[edge.dst]].first], [points[rev[edge.src]].second, points[rev[edge.dst]].second], linewidth=2, color="grey60")#, linealpha=0.5)
                    # display(Plots.plot!(title="Edge ($(rev[edge.src]), $(rev[edge.dst]))"))
                end
            end
            # display(Plots.plot!(title="Conflict Graph - Biclique $j"))
        end

        for edge in biclique_edges
            Plots.plot!([points[rev[edge.first]].first, points[rev[edge.second]].first], [points[rev[edge.first]].second, points[rev[edge.second]].second], linewidth=2, color=colors[(j % 3) + 1])
            # display(Plots.plot!(title="Edge ($(rev[edge.first]), $(rev[edge.second]))"))
        end

        # See plot_points
        x = map(point -> point.first, points)
        y = map(point -> point.second, points)
        for i = 1:length(points)
            if i in cover_vec[j].first
                Plots.scatter!([x[i]],[y[i]], color="red", markersize = 7) #, series_annotations=([Plots.text(string(i), :right, 8, "courier")]))
            elseif i in cover_vec[j].second
                Plots.scatter!([x[i]],[y[i]], color="blue", markersize = 7) #, series_annotations=([Plots.text(string(i), :right, 8, "courier")]))
            else
                Plots.scatter!([x[i]],[y[i]], color="grey35", markersize = 7) #, series_annotations=([Plots.text(string(i), :right, 8, "courier")]))
            end
        end

        display(Plots.plot!(title="Biclique $j"))
        savefig("$name $j.pdf")
    end
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

Plots and labels points given. Optional argument to plot a subset of vertices.
"""
function plot_points(points::Vector{Pair{Rational{Int64},Rational{Int64}}}; vertices::Dict{Int,Int}=Dict{Int,Int}(), with_labels=true)
    
    Plots.plot!(legend=false)

    if !isempty(vertices)
        points_filtered = []
        for v in keys(vertices)
            push!(points_filtered, points[v])
        end
        x = map(point -> point.first, points_filtered)
        y = map(point -> point.second, points_filtered)
        if with_labels
            Plots.scatter!(x,y, color="red", markersize=7, series_annotations=([Plots.text(string(x), :right, 10, "courier") for x in keys(vertices)]))
        else
            Plots.scatter!(x,y, color="red", markersize=7)
        end
    else
        x = map(point -> point.first, points)
        y = map(point -> point.second, points)
        if with_labels
            Plots.scatter!(x,y, color="red", markersize=7, series_annotations=([Plots.text(string(x), :right, 10, "courier") for x in 1:length(points)]))
        else
            Plots.scatter!(x,y, color="red", markersize=7)
        end
    end

end

"""
    plot_steps(obstacles, x, y, θ)
Description.
"""
function plot_steps(obstacles, x, y, θ)
    plot()
    ClutteredEnvPathOpt.plot_field(obstacles);
    scatter!(x[1:2:end], y[1:2:end], color="red", markersize=5, series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:2:length(x)]));
    scatter!(x[2:2:end], y[2:2:end], color="blue", markersize=5, series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 2:2:length(x)]));
    quiver!(x, y, quiver=(0.075 * cos.(θ), 0.075 * sin.(θ)))
    # display(plot!(title="Footsteps"))
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
