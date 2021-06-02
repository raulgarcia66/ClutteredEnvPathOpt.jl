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
    acc = [0] # accumulator

    _o, points, _g, faces = ClutteredEnvPathOpt.construct_graph(obstacles)

    for face in faces
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        polygon = Polyhedra.polyhedron(
            v,
            Polyhedra.DefaultLibrary{Float64}(Gurobi.Optimizer)
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

function get_M_A_b(obstacles)
    Ms = []
    As = []
    bs = []
    acc = [0]

    _o, points, _g, faces = ClutteredEnvPathOpt.construct_graph(obstacles)

    # Solver over the unitcell for now
    u = Polyhedra.convexhull([0,0],[0,1],[1,0],[1,1])
    unitcell = Polyhedra.polyhedron(u)

    for face in faces
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        polygon = Polyhedra.polyhedron(
            v,
            Polyhedra.DefaultLibrary{Float64}(Gurobi.Optimizer)
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
            # ?? #JuMP.@constraint(sub_model, sub_x in polygon)  # TODO "no method matching copy" error WHY DOES THIS ONLY WORK THE FIRST TIME
            JuMP.@objective(sub_model, Max, A * sub_x)

            stat = JuMP.optimize!(sub_model)

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
function solve_deits(obstacles, N, f1, f2, g, Q_g, Q_r, q_t, L, delta_f_max)
    # model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.20))
    model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer))

    
    # model has scalar variables x, y, theta, bin var t and d (2x1 decision), p (2x1 decision)
    JuMP.@variable(model, x[1:N])
    JuMP.@variable(model, y[1:N])
    JuMP.@variable(model, theta[1:N])

    JuMP.@variable(model, t[1:N], Bin)

    # Objective
    f = vec([[x[j], y[j], theta[j]] for j in 1:N])
    f_no_theta = vec([[x[j], y[j]] for j in 1:N])
    JuMP.@objective(
        model,
        Min,
        ((f[N] - g)' * Q_g * (f[N] - g)) + sum(q_t * t) + sum([(f[j + 1] - f[j])' * Q_r * (f[j + 1] - f[j]) for j in 1:(N-1)])
    )
    # JuMP.@objective(
    #     model,
    #     Min,
    #     ((f[N] - g)' * Q_g * (f[N] - g)) + sum(q_t * t)
    # )

    # Big M constraints
    #M, A, b, acc = get_M_A_b_easy(obstacles)
    M, A, b, acc = get_M_A_b(obstacles)

    num_safe_regions = length(acc)-1
    JuMP.@variable(model, z[1:N, 1:num_safe_regions], Bin)

    for j in 1:N
        for r in 1:(num_safe_regions)
            ids = (acc[r]+1):(acc[r+1])
            for i in ids
                # JuMP.@constraint(model, A[r, :]' * [x[j], y[j]] <= b[r] * z[j, r] + M[r] * (1 - z[j, r]))
                # JuMP.@constraint(model, A[r, :]' * [x[j], y[j]] <= b[r] * z[j, r] + 1 * (1 - z[j, r]))
                # println("$i")
                JuMP.@constraint(model, A[i, :]' * [x[j], y[j]] <= b[i] * z[j, r] + M[i] * (1 - z[j, r]))
            end
        end
        JuMP.@constraint(model, sum(z[j, :]) == 1)
    end

    # Reachability
    s = [piecewiselinear(model, theta[j], range(0, stop=(2 * pi), length=L), sin) for j in 1:N]
    c = [piecewiselinear(model, theta[j], range(0, stop=(2 * pi), length=L), cos) for j in 1:N]

    d1 = 0.5
    d2 = 0.5
    p1 = [0, 0.1]   
    p2 = [0, -0.1]
    for j in 2:N
        JuMP.@constraint(
            model,
            [
                d1,
                x[j] - x[j - 1] - c[j] * p1[1] + s[j] * p1[2],
                y[j] - y[j - 1] - s[j] * p1[1] - c[j] * p1[2]
            ] in SecondOrderCone()
        )

        JuMP.@constraint(
            model,
            [
                d2,
                x[j] - x[j - 1] - c[j] * p2[1] + s[j] * p2[2],
                y[j] - y[j - 1] - s[j] * p2[1] - c[j] * p2[2]
            ] in SecondOrderCone()
        )
    end

    # Set T to punish extra steps
    for j in 1:N
        if j % 2 == 1
            JuMP.@constraint(model, t[j] => {x[j] == f1[1]})
            JuMP.@constraint(model, t[j] => {y[j] == f1[2]})
            JuMP.@constraint(model, t[j] => {theta[j] == f1[3]})
        else
            JuMP.@constraint(model, t[j] => {x[j] == f2[1]})
            JuMP.@constraint(model, t[j] => {y[j] == f2[2]})
            JuMP.@constraint(model, t[j] => {theta[j] == f2[3]})
        end
    end

    # Max step distance
    for j in 2:N
        JuMP.@constraint(model, [delta_f_max; f_no_theta[j] - f_no_theta[j - 1]] in SecondOrderCone())
        # JuMP.@constraint(model, norm(f_no_theta[j] - f_no_theta[j - 1]) <= delta_f_max)
    end

    JuMP.@constraint(model, f[1] .== f1)
    JuMP.@constraint(model, f[2] .== f2)

    # Solve
    stat = JuMP.optimize!(model)

    return value.(x), value.(y), value.(theta)
end

function plot_steps(obstacles, x, y, theta)
    ClutteredEnvPathOpt.plot_field(obstacles)
    ClutteredEnvPathOpt.plot_lines(obstacles)
    ClutteredEnvPathOpt.plot_intersections(obstacles)
end