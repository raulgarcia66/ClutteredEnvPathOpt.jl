using ClutteredEnvPathOpt
using LightGraphs
import Polyhedra
# import GLPK
using Pipe
using Plots
using JuMP, Gurobi, PiecewiseLinearOpt

function get_M_A_b_easy(obstacles)
    Ms = []
    As = []
    bs = []
    acc = [0] # accumulator # try using reduce() instead

    _, points, _, _, free_faces = ClutteredEnvPathOpt.construct_graph(obstacles)

    for face in free_faces
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        polygon = Polyhedra.polyhedron(
            v,
            Polyhedra.DefaultLibrary{Rational{Int64}}(Gurobi.Optimizer)
        )
        halfspaces = Polyhedra.hrep(polygon).halfspaces
        temp = acc[end]
        push!(acc, temp + length(halfspaces))
        
        # for i in v.points
        #     println("$i $face")
        # end

        for hs in halfspaces
            A = hs.a'
            b = hs.β

            M = maximum(map(x -> A * x, v.points))

            push!(Ms, M)
            push!(As, A)
            push!(bs, b)
        end
    end
    # acc = acc[2:end]
    
    return vcat(Ms...), vcat(As...), vcat(bs...), acc
end

function get_M_A_b(points, free_faces)
    Ms = []
    As = []
    bs = []
    acc = [0] # try using reduce() instead

    # _, points, _, _, free_faces = ClutteredEnvPathOpt.construct_graph(obstacles)

    # Solve over the unitcell
    u = Polyhedra.convexhull([0//1,0//1],[0//1,1//1],[1//1,0//1],[1//1,1//1])
    unitcell = Polyhedra.polyhedron(u)

    for face in free_faces
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        polygon = Polyhedra.polyhedron(
            v,
            Polyhedra.DefaultLibrary{Rational{Int64}}(Gurobi.Optimizer)
        )
        halfspaces = Polyhedra.hrep(polygon).halfspaces
        temp = acc[end]
        push!(acc, temp + length(halfspaces))

        for hs in halfspaces
            A = hs.a'
            b = hs.β

            sub_model = JuMP.Model(Gurobi.Optimizer)
            #sub_x = JuMP.@variable(sub_model, [1:2])
            JuMP.@variable(sub_model, sub_x[1:2])
            JuMP.@constraint(sub_model, sub_x in unitcell)
            #JuMP.@constraint(sub_model, sub_x in safe_regions)
            JuMP.@objective(sub_model, Max, A * sub_x)

            JuMP.optimize!(sub_model)

            #return value.(sub_x)
            M = JuMP.objective_value(sub_model)

            push!(Ms, M)
            push!(As, A)
            push!(bs, b)
        end

        # A = vcat([hs.a for hs in halfspaces]...)
        # b = vcat([hs.β for hs in halfspaces]...)

        # M = map(hs ->
        #     begin
        #         A = hs.a
        #         b = hs.β

        #         sub_model = JuMP.Model(Gurobi.Optimizer)
        #         #sub_x = JuMP.@variable(sub_model, [1:2])
        #         JuMP.@variable(sub_model, sub_x[1:2])
        #         JuMP.@constraint(sub_model, sub_x in polygon)   # TODO "no method matching copy" error WHY DOES THIS ONLY WORK THE FIRST TIME
        #         JuMP.@objective(sub_model, Max, A' * sub_x)

        #         stat = JuMP.optimize!(sub_model)

        #         #return value.(sub_x)
        #         return JuMP.objective_value(sub_model)
        #     end,
        #     halfspaces
        # )

        # push!(Ms, M)
        # push!(As, A)
        # push!(bs, b)
    end
    # acc = acc[2:end]

    return vcat(Ms...), vcat(As...), vcat(bs...), acc
    #return Ms, As, bs, acc
end

function plot_steps(obstacles, x, y, θ)
    plot()
    ClutteredEnvPathOpt.plot_field(obstacles);
    scatter!(x[1:2:end], y[1:2:end], color="red", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:2:length(x)]));
    scatter!(x[2:2:end], y[2:2:end], color="blue", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 2:2:length(x)]));
    quiver!(x, y, quiver=(0.075 * cos.(θ), 0.075 * sin.(θ)))
    # display(plot!(title="Footsteps"))
end

function plot_circles(x, y, theta; R1=0.20, R2=0.20, p1=[0, 0.07], p2=[0, -0.27])
    angles = LinRange(0,2*pi, 100)
    #plot()
    for i = 2:length(x)
        plot()
        if i % 2 == 1
            circlex1 = R1*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(-p1[1]) - sin(theta[i-1])*(-p1[2]) )
            circley1 = R1*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(-p1[1]) + cos(theta[i-1])*(-p1[2]) )

            circlex2 = R2*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(-p2[1]) - sin(theta[i-1])*(-p2[2]) )
            circley2 = R2*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(-p2[1]) + cos(theta[i-1])*(-p2[2]) )
            #scatter!([x[i-1]; x[i]], [y[i-1]; y[i]])
            scatter!([x[i-1]], [y[i-1]], color = "blue")
            scatter!([x[i]], [y[i]], color = "red")
            quiver!([x[i-1]; x[i]], [y[i-1]; y[i]], 
                quiver=([0.075 * cos(theta[i-1]); 0.075 * cos(theta[i])], 
                [0.075 * sin(theta[i-1]); 0.075 * sin(theta[i])])
                )
        else
            circlex1 = R1*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(p1[1]) - sin(theta[i-1])*(p1[2]) )
            circley1 = R1*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(p1[1]) + cos(theta[i-1])*(p1[2]) )

            circlex2 = R2*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(p2[1]) - sin(theta[i-1])*(p2[2]) )
            circley2 = R2*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(p2[1]) + cos(theta[i-1])*(p2[2]) )
            #scatter!([x[i-1]; x[i]], [y[i-1]; y[i]])
            scatter!([x[i-1]], [y[i-1]], color = "red")
            scatter!([x[i]], [y[i]], color = "blue")
            quiver!([x[i-1]; x[i]], [y[i-1]; y[i]], 
                quiver=([0.075 * cos(theta[i-1]); 0.075 * cos(theta[i])], 
                [0.075 * sin(theta[i-1]); 0.075 * sin(theta[i])])
                )
        end
        plot!(circlex1, circley1, color="dodgerblue")
        plot!(circlex2, circley2, color="maroon")
        display(plot!(legend=false, xlims=(-0.1,1.1), ylims=(-0.1,1.1)))
    end
    # display(plot!(legend=false, xlims=(-0.1,1.1), ylims=(-0.1,1.1)))
end

# obstacles <- obstacles (list of polyhedra)
# N <- max number of steps (scalar)
# f1 <- initial left foot ([x, y, theta])
# f2 <- initial right foot ([x, y, theta])
# g <- goal pose
# Q_g <- pose weights to goal (4x4 mat)
# Q_r <- pose weights between steps (4x4 mat)
# q_t <- weight on unused steps (scalar)
# method <- "merged" for the compact biclique cover, "full" for the original biclique cover, "bigM" for big-M constraints
# d1 = 0.2 <- radius of reference foot circle
# d2 = 0.2 <- radius of moving foot circle
# p1 = [0, 0.07] <- center of reference foot circle
# p2 = [0, -0.27] <- center of moving foot circle
# delta_x_y_max = 0.1 <- max stride norm in space
# delta_θ_max = pi/4 <- max difference in θ
# L = 5 <- number of pieces of pwl sin/cos (scalar)
"""
    solve_steps()

Adapted from Deits and Tedrake 2014.
"""
function solve_steps(obstacles, N, f1, f2, g, Q_g, Q_r, q_t; method="merged", partition="CDT", merge_faces=true, d1=0.2, d2=0.2, p1=[0, 0.07], p2=[0, -0.27], relax=false)

    if partition == "CDT"
        _, points, graph, _, free_faces = ClutteredEnvPathOpt.construct_graph_delaunay(obstacles, merge_faces=merge_faces)
    else
        _, points, graph, _, free_faces = ClutteredEnvPathOpt.construct_graph(obstacles)
    end

    skeleton = LabeledGraph(graph)

    # model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer))
    # model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer, "Heuristics"=> 0, "Cuts"=> 0, "Precrush"=>1, "MIPGap" => .01, "TimeLimit" => 300))
    if relax
        model = JuMP.Model(Gurobi.Optimizer)
    else
        model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => .01, "TimeLimit" => 300))
    end

    # model has scalar variables x, y, θ, binary variable t
    JuMP.@variable(model, x[1:N])
    JuMP.@variable(model, y[1:N])
    JuMP.@variable(model, θ[1:N])
    JuMP.@variable(model, t[1:N], Bin)

    # Objective
    f = vec([[x[j], y[j], θ[j]] for j in 1:N])
    # f_no_theta = vec([[x[j], y[j]] for j in 1:N])
    JuMP.@objective(
        model,
        Min,
        ((f[N] - g)' * Q_g * (f[N] - g)) + sum(q_t * t) + sum((f[j + 1] - f[j])' * Q_r * (f[j + 1] - f[j]) for j in 1:(N-1))
    )

    cover = Set{Pair{Set{Int64}, Set{Int64}}}()
    merged_cover = Set{Pair{Set{Int64}, Set{Int64}}}()
    num_free_face_ineq = 0

    if method != "bigM"
        cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces)
        feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
        merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)

        (valid_cover, _, _, _) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover)
        (valid_cover_merged, _, _, _) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover)
    end
    
    # Footstep location constraints
    if method == "merged" && valid_cover_merged
        J = LightGraphs.nv(skeleton.graph)
        JuMP.@variable(model, λ[1:N, 1:J] >= 0)
        JuMP.@variable(model, z[1:N, 1:length(merged_cover)], Bin)
        for i in 1:N
            for (j,(A,B)) in enumerate(merged_cover)
                JuMP.@constraint(model, sum(λ[i, v] for v in A) <= z[i, j])
                JuMP.@constraint(model, sum(λ[i, v] for v in B) <= 1 - z[i, j])
            end
            JuMP.@constraint(model, sum(λ[i,:]) == 1)
            JuMP.@constraint(model, x[i] == sum(points[j].first * λ[i, j] for j in 1:J))
            JuMP.@constraint(model, y[i] == sum(points[j].second * λ[i, j] for j in 1:J))
        end
    elseif method != "bigM" && valid_cover
        J = LightGraphs.nv(skeleton.graph)
        JuMP.@variable(model, λ[1:N, 1:J] >= 0)
        JuMP.@variable(model, z[1:N, 1:length(cover)], Bin)
        for i in 1:N
            for (j,(A,B)) in enumerate(cover)
                JuMP.@constraint(model, sum(λ[i, v] for v in A) <= z[i, j])
                JuMP.@constraint(model, sum(λ[i, v] for v in B) <= 1 - z[i, j])
            end
            JuMP.@constraint(model, sum(λ[i,:]) == 1)
            JuMP.@constraint(model, x[i] == sum(points[j].first * λ[i, j] for j in 1:J))
            JuMP.@constraint(model, y[i] == sum(points[j].second * λ[i, j] for j in 1:J))
        end
    else
        # Runs if method = "bigM", or if a valid biclique cover is not found
        M, A, b, acc = get_M_A_b(points, free_faces)
        num_free_face_ineq = length(M)
        JuMP.@variable(model, z[1:N, 1:length(free_faces)], Bin)
        for j in 1:N
            for r in 1:length(free_faces)
                ids = (acc[r]+1):(acc[r+1])
                for i in ids
                    # JuMP.@constraint(model, A[r, :]' * [x[j], y[j]] <= b[r] * z[j, r] + M[r] * (1 - z[j, r]))
                    # JuMP.@constraint(model, A[r, :]' * [x[j], y[j]] <= b[r] * z[j, r] + 1 * (1 - z[j, r]))
                    # println("$i")
                    JuMP.@constraint(model, A[i, :]' * [x[j]; y[j]] <= b[i] * z[j, r] + M[i] * (1 - z[j, r]))
                end
            end
            JuMP.@constraint(model, sum(z[j, :]) == 1)
        end
    end

    # Reachability
    # Breakpoints need to be strategically chosen
    s_break_pts = [0, 5pi/16, 11pi/16, 21pi/16, 27pi/16, 2pi]
    c_break_pts = [0, 3pi/16, 13pi/16, 19pi/16, 29pi/16, 2pi]
    s = [piecewiselinear(model, θ[j], s_break_pts, sin) for j in 1:N]
    c = [piecewiselinear(model, θ[j], c_break_pts, cos) for j in 1:N]

    # For footstep j, the cirles are created from the frame of reference of footstep j-1
    # Use negative if in frame of reference of right foot (j = 1 is left foot)
    # Position is with respect to θ = 0
    for j in 2:N
        if j % 2 == 1
            # In frame of reference of right foot
            JuMP.@constraint(
                model,
                [
                    d1,
                    x[j] - x[j - 1] - c[j-1] * (-p1[1]) + s[j-1] * (-p1[2]),
                    y[j] - y[j - 1] - s[j-1] * (-p1[1]) - c[j-1] * (-p1[2])
                ] in SecondOrderCone()
            )

            JuMP.@constraint(
                model,
                [
                    d2,
                    x[j] - x[j - 1] - c[j-1] * (-p2[1]) + s[j-1] * (-p2[2]),
                    y[j] - y[j - 1] - s[j-1] * (-p2[1]) - c[j-1] * (-p2[2])
                ] in SecondOrderCone()
            )
        else
            # In frame of reference left foot
            JuMP.@constraint(
                model,
                [
                    d1,
                    x[j] - x[j - 1] - c[j-1] * p1[1] + s[j-1] * p1[2],
                    y[j] - y[j - 1] - s[j-1] * p1[1] - c[j-1] * p1[2]
                ] in SecondOrderCone()
            )

            JuMP.@constraint(
                model,
                [
                    d2,
                    x[j] - x[j - 1] - c[j-1] * p2[1] + s[j-1] * p2[2],
                    y[j] - y[j - 1] - s[j-1] * p2[1] - c[j-1] * p2[2]
                ] in SecondOrderCone()
            )
        end
    end

    # Set T to punish extra steps
    for j in 1:N
        if j % 2 == 1
            JuMP.@constraint(model, t[j] => {x[j] == f1[1]}) # use f1 and f2?
            JuMP.@constraint(model, t[j] => {y[j] == f1[2]})
            JuMP.@constraint(model, t[j] => {θ[j] == f1[3]})
        else
            JuMP.@constraint(model, t[j] => {x[j] == f2[1]})
            JuMP.@constraint(model, t[j] => {y[j] == f2[2]})
            JuMP.@constraint(model, t[j] => {θ[j] == f2[3]})
        end
    end

    # Max step distance
    # for j in 3:2:(N-1)
    #     JuMP.@constraint(model, [delta_x_y_max; f_no_theta[j] - f_no_theta[j - 2]] in SecondOrderCone())
    #     JuMP.@constraint(model, [delta_x_y_max; f_no_theta[j + 1] - f_no_theta[j - 1]] in SecondOrderCone())
    # end

    # Max theta difference for individual feet
    # for j in 3:2:(N-1)
    #     JuMP.@constraint(model, -delta_θ_max <= θ[j] - θ[j - 2] <= delta_θ_max)
    #     JuMP.@constraint(model, -delta_θ_max <= θ[j + 1] - θ[j - 1] <= delta_θ_max)
    # end

    # Max θ difference between feet
    for j in 2:N
        JuMP.@constraint(model, -pi/8 <= θ[j] - θ[j-1] <= pi/8)
    end

    # Initial footstep positions
    JuMP.@constraint(model, f[1] .== f1)
    JuMP.@constraint(model, f[2] .== f2)

    # Solve
    if relax
        relax_integrality(model)
    end
    JuMP.optimize!(model)

    if method == "merged" && valid_cover_merged
        println("\n\nUsed merged cover.\n\n")
    elseif method != "bigM" && valid_cover 
        println("\n\nUsed full cover.\n\n")
        method = "full"
    else
        println("\n\nUsed big-M constraints.\n\n")
        method = "bigM"
    end

    # TODO: CHECK IF MODEL HAS VALUES FIRST?
    stats = (termination_status(model), objective_value(model), solve_time(model), relative_gap(model), simplex_iterations(model), 
            node_count(model), LightGraphs.nv(skeleton.graph), length(merged_cover), length(cover), length(free_faces), num_free_face_ineq, method)

    # TODO: SAVE PLOTS OF SOLUTIONS OF MIQCQP
    # MOVE THIS OUTSIDE FUNCTION?
    # if !relax && has_values(model)
    #     x = value.(x)
    #     y = value.(y)
    #     θ = value.(θ)
    #     t = value.(t)
    #     num_to_trim = length(filter(tj -> tj > 0.5, t[3:end]))
    #     if num_to_trim % 2 == 0
    #         x = vcat(x[1:2], x[num_to_trim + 3 : end]);
    #         y = vcat(y[1:2], y[num_to_trim + 3 : end]);
    #         θ = vcat(θ[1:2], θ[num_to_trim + 3 : end]);
    #     else
    #         x = vcat(x[1], x[num_to_trim + 3 : end]);
    #         y = vcat(y[1], y[num_to_trim + 3 : end]);
    #         θ = vcat(θ[1], θ[num_to_trim + 3 : end]);
    #     end
    #     plot_steps(obstacles, x, y, θ)
    #     png("IMG_NAME")
    # end

    return value.(x), value.(y), value.(θ), value.(t), stats
end
