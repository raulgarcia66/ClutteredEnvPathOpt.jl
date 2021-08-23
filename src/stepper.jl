using ClutteredEnvPathOpt
using LightGraphs
import Polyhedra
import GLPK
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

            M = maximum(map(x -> A * x, v.points)) # think it needs to be over all subfaces

            push!(Ms, M)
            push!(As, A)
            push!(bs, b)
        end
    end
    # acc = acc[2:end]
    
    return vcat(Ms...), vcat(As...), vcat(bs...), acc
end

function get_M_A_b(obstacles)
    Ms = []
    As = []
    bs = []
    acc = [0] # try using reduce() instead

    _, points, _, _, free_faces = ClutteredEnvPathOpt.construct_graph(obstacles)

    # Solve over the unitcell for now
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

# obstacles <- obstacles (list of polyhedra, see email)
# N <- max number of steps (scalar)
# f1 <- initial left foot ([x, y, theta])
# f2 <- initial right foot ([x, y, theta])
# g <- goal pose
# Q_g <- pose weights to goal (4x4 mat) (just identity for now)
# Q_r <- pose weights between steps (4x4 mat)
# q_t <- weight on unused steps (scalar)
# L <- number of pieces of pwl sin/cos (scalar)
# delta_f_max <- max stride norm
function solve_deits(obstacles, N, f1, f2, g, Q_g, Q_r, q_t, L, delta_f_max; d1=0.2, d2=0.2, p1=[0, 0.05], p2=[0, -0.25])

    _, points, graph, obstacle_faces, free_faces = ClutteredEnvPathOpt.construct_graph(obstacles)
    skeleton = LabeledGraph(graph)
    # all_faces = union(obstacle_faces, free_faces)

    model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => .01))
    # "TimeLimit" => 150, "MIPGap" => 0.20
    
    # model has scalar variables x, y, theta, bin var t and d (2x1 decision), p (2x1 decision)
    JuMP.@variable(model, x[1:N])
    JuMP.@variable(model, y[1:N])
    JuMP.@variable(model, θ[1:N])

    JuMP.@variable(model, t[1:N], Bin)

    # Objective
    f = vec([[x[j], y[j], θ[j]] for j in 1:N])
    f_no_theta = vec([[x[j], y[j]] for j in 1:N])
    JuMP.@objective(
        model,
        Min,
        ((f[N] - g)' * Q_g * (f[N] - g)) + sum(q_t * t) + sum((f[j + 1] - f[j])' * Q_r * (f[j + 1] - f[j]) for j in 1:(N-1))
    )

    # # Big M constraints for footstep location
    # M, A, b, acc = get_M_A_b(obstacles)
    # num_free_faces = length(acc)-1
    # JuMP.@variable(model, z[1:N, 1:num_free_faces], Bin)
    # for j in 1:N
    #     for r in 1:num_free_faces
    #         ids = (acc[r]+1):(acc[r+1])
    #         for i in ids
    #             # JuMP.@constraint(model, A[r, :]' * [x[j], y[j]] <= b[r] * z[j, r] + M[r] * (1 - z[j, r]))
    #             # JuMP.@constraint(model, A[r, :]' * [x[j], y[j]] <= b[r] * z[j, r] + 1 * (1 - z[j, r]))
    #             # println("$i")
    #             JuMP.@constraint(model, A[i, :]' * [x[j]; y[j]] <= b[i] * z[j, r] + M[i] * (1 - z[j, r]))
    #         end
    #     end
    #     JuMP.@constraint(model, sum(z[j, :]) == 1)
    # end

    # Disjunctive constraints for footstep location
    cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces)
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
    (valid_cover, _, missing_edges, _) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover)
    
    if valid_cover
        J = LightGraphs.nv(skeleton.graph)
        JuMP.@variable(model, λ[1:N, 1:J] >= 0)
        JuMP.@variable(model, z[1:N, 1:length(cover)], Bin)
        for i = 1:N
            for (j,(A,B)) in enumerate(cover)
                JuMP.@constraint(model, sum(λ[i, v] for v in A) <= z[i, j])
                JuMP.@constraint(model, sum(λ[i, v] for v in B) <= 1 - z[i, j])
            end
            JuMP.@constraint(model, sum(λ[i,:]) == 1)
            JuMP.@constraint(model, x[i] == sum(points[j].first * λ[i, j] for j in 1:J))
            JuMP.@constraint(model, y[i] == sum(points[j].second * λ[i, j] for j in 1:J))
        end
    else
        error("Valid Cover Not Found")
        # # Big M constraints
        # M, A, b, acc = get_M_A_b(obstacles)
        # num_free_faces = length(acc)-1
        # JuMP.@variable(model, z[1:N, 1:num_free_faces], Bin)
        # for j in 1:N
        #     for r in 1:num_free_faces
        #         ids = (acc[r]+1):(acc[r+1])
        #         for i in ids
        #             # JuMP.@constraint(model, A[r, :]' * [x[j], y[j]] <= b[r] * z[j, r] + M[r] * (1 - z[j, r]))
        #             # JuMP.@constraint(model, A[r, :]' * [x[j], y[j]] <= b[r] * z[j, r] + 1 * (1 - z[j, r]))
        #             # println("$i")
        #             JuMP.@constraint(model, A[i, :]' * [x[j]; y[j]] <= b[i] * z[j, r] + M[i] * (1 - z[j, r]))
        #         end
        #     end
        #     JuMP.@constraint(model, sum(z[j, :]) == 1)
        # end
    end

    # # Alternately with anonymous variabels
    # if valid_cover
    #     J = LightGraphs.nv(skeleton.graph)
    #     for i = 1:N
    #         λ = JuMP.@variable(model, [1:J], lower_bound = 0)
    #         z = JuMP.@variable(model, [1:length(cover)], Bin)
    #         for (j,(A,B)) in enumerate(cover)
    #             JuMP.@constraint(model, sum(λ[v] for v in A) <= z[j])
    #             JuMP.@constraint(model, sum(λ[v] for v in B) <= 1 - z[j])
    #         end
    #         JuMP.@constraint(model, sum(λ) == 1)
    #         JuMP.@constraint(model, x[i] == sum(points[j].first * λ[j] for j in 1:J))
    #         JuMP.@constraint(model, y[i] == sum(points[j].second * λ[j] for j in 1:J))
    #     end
    # else
    #     error("Valid Cover Not Found")
    #     # # Big M constraints
    #     # M, A, b, acc = get_M_A_b(obstacles)
    #     # num_free_faces = length(acc)-1
    #     # JuMP.@variable(model, z[1:N, 1:num_free_faces], Bin)
    #     # for j in 1:N
    #     #     for r in 1:num_free_faces
    #     #         ids = (acc[r]+1):(acc[r+1])
    #     #         for i in ids
    #     #             # JuMP.@constraint(model, A[r, :]' * [x[j], y[j]] <= b[r] * z[j, r] + M[r] * (1 - z[j, r]))
    #     #             # JuMP.@constraint(model, A[r, :]' * [x[j], y[j]] <= b[r] * z[j, r] + 1 * (1 - z[j, r]))
    #     #             # println("$i")
    #     #             JuMP.@constraint(model, A[i, :]' * [x[j]; y[j]] <= b[i] * z[j, r] + M[i] * (1 - z[j, r]))
    #     #         end
    #     #     end
    #     #     JuMP.@constraint(model, sum(z[j, :]) == 1)
    #     # end
    # end

    # Reachability
    s = [piecewiselinear(model, θ[j], range(0, stop=(2 * pi), length=L), sin) for j in 1:N]
    c = [piecewiselinear(model, θ[j], range(0, stop=(2 * pi), length=L), cos) for j in 1:N]

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
            JuMP.@constraint(model, t[j] => {x[j] == f1[1]})
            JuMP.@constraint(model, t[j] => {y[j] == f1[2]})
            JuMP.@constraint(model, t[j] => {θ[j] == f1[3]})
        else
            JuMP.@constraint(model, t[j] => {x[j] == f2[1]})
            JuMP.@constraint(model, t[j] => {y[j] == f2[2]})
            JuMP.@constraint(model, t[j] => {θ[j] == f2[3]})
        end
    end

    # Max step distance
    for j in 2:N
        JuMP.@constraint(model, [delta_f_max; f_no_theta[j] - f_no_theta[j - 1]] in SecondOrderCone())
    end

    # Initial footstep positions
    JuMP.@constraint(model, f[1] .== f1)
    JuMP.@constraint(model, f[2] .== f2)

    # Solve
    JuMP.optimize!(model)

    return value.(x), value.(y), value.(θ), value.(t)
end

function plot_steps(obstacles, x, y, theta)
    plot()
    ClutteredEnvPathOpt.plot_field(obstacles);
    scatter!(x[1:2:end], y[1:2:end], color="red", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 1:2:length(x)]));
    scatter!(x[2:2:end], y[2:2:end], color="blue", series_annotations=([Plots.text(string(x), :right, 8, "courier") for x in 2:2:length(x)]));
    quiver!(x, y, quiver=(0.075 * cos.(theta), 0.075 * sin.(theta)))
    display(plot!(title="Footsteps"))
end

function plot_circles(R1, R2, p1, p2, x, y, theta)
    angles = LinRange(0,2*pi, 100)
    #plot()
    for i = 2:length(x)
        plot()
        if i % 2 == 1
            circlex1 = R1*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(-p1[1]) - sin(theta[i-1])*(-p1[2]) )
            circley1 = R1*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(-p1[1]) + cos(theta[i-1])*(-p1[2]) )

            circlex2 = R2*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(-p2[1]) - sin(theta[i-1])*(-p2[2]) )
            circley2 = R2*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(-p2[1]) + cos(theta[i-1])*(-p2[2]) )
            scatter!([x[i-1]; x[i]], [y[i-1]; y[i]])
        else
            circlex1 = R1*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(p1[1]) - sin(theta[i-1])*(p1[2]) )
            circley1 = R1*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(p1[1]) + cos(theta[i-1])*(p1[2]) )

            circlex2 = R2*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(p2[1]) - sin(theta[i-1])*(p2[2]) )
            circley2 = R2*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(p2[1]) + cos(theta[i-1])*(p2[2]) )
            scatter!([x[i-1]; x[i]], [y[i-1]; y[i]])
        end
        plot!(circlex1, circley1, color="dodgerblue")
        plot!(circlex2, circley2, color="maroon")
        display(plot!(legend=false, xlims=(-0.1,1.1), ylims=(-0.1,1.1)))
    end
    # display(plot!(legend=false, xlims=(-0.1,1.1), ylims=(-0.1,1.1)))
end