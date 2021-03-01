using ClutteredEnvPathOpt
using LightGraphs
using Pipe
using Plots
using JuMP, Gurobi

function get_M(obstacles)
    for obstacle in obstacles
        halfspaces = obstacle.hrep.halfspaces

        
    end
end

# M <- obstacles (matrix, see email)
# N <- max number of steps (scalar)
# f1 <- initial left foot ([x, y, z, theta])
# f2 <- initial right foot ([x, y, z, theta])
# g <- goal pose
# Q_g <- pose weights to goal (4x4 mat) (just identity for now)
# Q_r <- pose weights between steps (4x4 mat)
# q_t <- weight on unused steps (scalar)
# L <- number of pieces of pwl sin/cos (scalar)
function solve_deits(M, N, f1, f2, g, Q_g, Q_r, q_t, L)
    model = JuMP.Model(Gurobi.Optimizer)
    
    # model has scalar variables x, y, z?, theta and d (2x1 decision), p (2x1 decision)
    JuMP.@variable(model, x[1:N])
    JuMP.@variable(model, y[1:N])
    JuMP.@variable(model, z[1:N])
    JuMP.@variable(model, theta[1:N])

    JuMP.@variable(model, d1[1:N])
    JuMP.@variable(model, d2[1:N])

    JuMP.@variable(model, p1[1:N, 1:2])
    JuMP.@variable(model, p2[1:N, 1:2])

    # Objective
    # TODO: How do I express j in a julia sum function
    JuMP.@objective(
        model,
        Min,
        (([x[N], y[N], z[N], theta[N]] - g)' * Q_g * ([x[N], y[N], z[N], theta[N]] - g)) + sum([]) + sum([])
    )

    # Big M
    # TODO 

    # Reachability
    s = [piecewiselinear(model, theta[j], range(0, stop=(2 * pi), length=L), sin) for j in 1:N]
    c = [piecewiselinear(model, theta[j], range(0, stop=(2 * pi), length=L), cos) for j in 1:N]

    # Start on second step?
    for j in 2:N
        JuMP.@constraint(
            model,
            [
                d1[j],
                x[j] - x[j - 1] - c[j] * p1[j, 1] + s[j] * p1[j, 2],
                y[j] - y[j - 1] - s[j] * p1[j, 1] - c[j] * p1[j, 2]
            ] in SecondOrderCone()
        )

        JuMP.@constraint(
            model,
            [
                d2[j],
                x[j] - x[j - 1] - c[j] * p2[j, 1] + s[j] * p2[j, 2],
                y[j] - y[j - 1] - s[j] * p2[j, 1] - c[j] * p2[j, 2]
            ] in SecondOrderCone()
        )
    end

    # Fix initial pose
    JuMP.@constraint(model, x[1] = f1[1])
    JuMP.@constraint(model, y[1] = f1[2])
    JuMP.@constraint(model, z[1] = f1[3])
    JuMP.@constraint(model, theta[1] = f1[4])

    JuMP.@constraint(model, x[2] = f2[1])
    JuMP.@constraint(model, y[2] = f2[2])
    JuMP.@constraint(model, z[2] = f2[3])
    JuMP.@constraint(model, theta[2] = f2[4])

    # Solve
end