using ClutteredEnvPathOpt
using LightGraphs
using Pipe
using Plots
using JuMP, Gurobi

function get_M_A_b(obstacles)
    Ms = []
    As = []
    bs = []

    for obstacle in obstacles
        halfspaces = obstacle.hrep.halfspaces
        A = hcat([hs.a for hs in halfspaces]...)'
        b = [hs.β for hs in halfspaces]

        push!(Ms, 1) # TODO
        push!(As, A)
        push!(bs, b)
    end

    return Ms, As, bs
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
function solve_deits(obstacles, N, f1, f2, g, Q_g, Q_r, q_t, L)
    model = JuMP.Model(Gurobi.Optimizer)
    
    # model has scalar variables x, y, z?, theta and d (2x1 decision), p (2x1 decision)
    JuMP.@variable(model, x[1:N])
    JuMP.@variable(model, y[1:N])
    JuMP.@variable(model, z[1:N])   # QUESTION: Do we need this variable since we're in the plane? We don't seem to use it anywhere
    JuMP.@variable(model, theta[1:N])

    JuMP.@variable(model, d1[1:N])
    JuMP.@variable(model, d2[1:N])

    JuMP.@variable(model, p1[1:N, 1:2])
    JuMP.@variable(model, p2[1:N, 1:2])

    # Objective
    # QUESTION: How do I express these sums here since I need to do some sort of indexing? Typically I would just call sum() but I dont think we can here
    JuMP.@objective(
        model,
        Min,
        (([x[N], y[N], z[N], theta[N]] - g)' * Q_g * ([x[N], y[N], z[N], theta[N]] - g)) + sum([]) + sum([])
    )

    # Big M constraints
    M, A, b = get_M_A_b(obstacles)
    JuMP.@variable(model, zzz[1:N, 1:length(obstacles)], Bin)

    for j in 1:N
        for r in 1:length(obstacles)
            JuMP.@constraint(model, A[r] * [x[j], y[j]] <= b[r] * zzz[j, r] + M[r] * (1 - zzz[j, r]))
        end

        JuMP.@constraint(model, sum(zzz[j, :]) = 1)
    end

    # Reachability
    s = [piecewiselinear(model, theta[j], range(0, stop=(2 * pi), length=L), sin) for j in 1:N]
    c = [piecewiselinear(model, theta[j], range(0, stop=(2 * pi), length=L), cos) for j in 1:N]

    # QUESTION: Is it correct to start on second step here?
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
    stat = JuMP.optimize!(model)

    return stat
end