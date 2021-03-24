using ClutteredEnvPathOpt
using LightGraphs
import Polyhedra
import GLPK
using Pipe
using Plots
using JuMP, Gurobi, PiecewiseLinearOpt

function get_M_A_b(obstacles)
    Ms = []
    As = []
    bs = []

    _o, points, _g, faces = ClutteredEnvPathOpt.construct_graph(obstacles)

    for face in faces
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)
        polygon = Polyhedra.polyhedron(
            v,
            Polyhedra.DefaultLibrary{Float64}(Gurobi.Optimizer)
        )
        halfspaces = Polyhedra.hrep(polygon).halfspaces

        A = hcat([hs.a for hs in halfspaces]...)'
        b = [hs.β for hs in halfspaces]

        M = map(hs ->
            begin
                A = hs.a
                b = hs.β

                sub_model = JuMP.Model(Gurobi.Optimizer)
                JuMP.@variable(sub_model, sub_x[1:2])
                JuMP.@constraint(sub_model, sub_x in polygon)   # TODO "no method matching copy" error WHY DOES THIS ONLY WORK THE FIRST TIME
                JuMP.@objective(sub_model, Max, A' * sub_x)

                stat = JuMP.optimize!(sub_model)

                return value.(sub_x)
            end,
            halfspaces
        )

        push!(Ms, M)
        push!(As, A)
        push!(bs, b)
    end

    return Ms, As, bs
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
    model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer))
    
    # model has scalar variables x, y, theta, bin var t and d (2x1 decision), p (2x1 decision)
    JuMP.@variable(model, x[1:N])
    JuMP.@variable(model, y[1:N])
    JuMP.@variable(model, theta[1:N])

    JuMP.@variable(model, t[1:N], Bin)

    # JuMP.@variable(model, d1[1:N])
    # JuMP.@variable(model, d2[1:N])

    # JuMP.@variable(model, p1[1:N, 1:2])
    # JuMP.@variable(model, p2[1:N, 1:2])

    # Objective
    f = vec([[x[j], y[j], theta[j]] for j in 1:N])
    JuMP.@objective(
        model,
        Min,
        ((f[N] - g)' * Q_g * (f[N] - g)) + sum(q_t * t) + sum([(f[j + 1] - f[j])' * Q_r * (f[j + 1] - f[j]) for j in 1:(N-1)])
    )

    # Big M constraints
    # M, A, b = get_M_A_b(obstacles)
    # JuMP.@variable(model, z[1:N, 1:length(obstacles)], Bin)

    # for j in 1:N
    #     for r in 1:length(obstacles)
    #         JuMP.@constraint(model, A[r] * [x[j], y[j]] <= b[r] * z[j, r] + M[r] * (1 - z[j, r]))
    #     end

    #     JuMP.@constraint(model, sum(z[j, :]) == 1)
    # end

    # Reachability
    s = [piecewiselinear(model, theta[j], range(0, stop=(2 * pi), length=L), sin) for j in 1:N]
    c = [piecewiselinear(model, theta[j], range(0, stop=(2 * pi), length=L), cos) for j in 1:N]

    d1 = 1
    d2 = 1
    p1 = 0
    p2 = 0
    for j in 2:N
        JuMP.@constraint(
            model,
            [
                1,
                x[j] - x[j - 1] - c[j] * p1 + s[j] * p1,
                y[j] - y[j - 1] - s[j] * p1 - c[j] * p1
            ] in SecondOrderCone()
        )

        JuMP.@constraint(
            model,
            [
                1,
                x[j] - x[j - 1] - c[j] * p2 + s[j] * p2,
                y[j] - y[j - 1] - s[j] * p2 - c[j] * p2
            ] in SecondOrderCone()
        )
    end

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

    # Max step distance TODO: did I do second order cones correct here?
    for j in 2:N
        # JuMP.@constraint(model, norm(f[j] - f[j - 1]) <= delta_f_max)
        # @constraint(model, norm(x) <= t)` should now be written as `@constraint(model, [t; x] in SecondOrderCone())
        JuMP.@constraint(model, [delta_f_max; f[j] - f[j - 1]] in SecondOrderCone())

    end

    # Solve
    stat = JuMP.optimize!(model)

    return value.(x), value.(y), value.(theta)
end

function plot_steps(obstacles, x, y, theta)
    ClutteredEnvPathOpt.plot_field(obstacles)
    ClutteredEnvPathOpt.plot_lines(obstacles)
    ClutteredEnvPathOpt.plot_intersections(obstacles)
end