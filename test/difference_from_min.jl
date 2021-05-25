#ENV["GUROBI_HOME"] = "/Library/gurobi910/mac64"

using ClutteredEnvPathOpt
using LightGraphs
using Pipe
using Plots
using JuMP, Gurobi

function plot_scatter(num_obstacles, num_samples)
    io = open("points.txt", "w")
    close(io)

    X = []
    Y = []
    N = []

    for i in 1:num_samples
        println("SAMPLE: ", i)

        x, y, n = gen_point(num_obstacles)

        push!(X, x)
        push!(Y, y)
        push!(N, n)
    end

    scatter(X, Y, title="Computed Biclique Cover Sizes vs. Optimal")
    xlabel!("Computed Biclique Cover Sizes")
    ylabel!("Optimal (via ILP)")
    Plots.savefig("scatter")

    io = open("points.txt", "a")
    println(io, "X=", X)
    println(io, "Y=", Y)
    println(io, "N=", N)
    close(io)

    return (X, Y, N)
end

function gen_point(num_obstacles)
    obstacles, points, g, faces = ClutteredEnvPathOpt.plot_new(num_obstacles)
    skeleton = LabeledGraph(g)

    cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, faces)
    x = length(cover)

    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(faces))
    E = @pipe complement(feg.graph) |> incidence_matrix(_)

    lower = Int(ceil(log2(length(faces))))
    upper = x

    while upper - lower > 1
        t = (upper + lower) ÷ 2 # integer division
        is_feasible = is_biclique_feasible(E, t)
        
        if is_feasible
            upper = t
        else
            lower = t
        end
    end

    y = lower
    plot_optimal(E, obstacles, points, lower, upper)

    io = open("points.txt", "a")
    println(io, "x=",x,";y=",y,";")
    close(io)

    return (x, y, length(points))
end

function is_biclique_feasible(E, t)
    println("CHECKING FEASIBILITY: t = ", t)

    n = size(E)[1]
    J = 1:n

    model = JuMP.Model(Gurobi.Optimizer)
    JuMP.@variable(model, x[1:t,1:n], Bin)
    JuMP.@variable(model, y[1:t,1:n], Bin)
    JuMP.@variable(model, z[1:t,1:n,1:n], Bin)

    for j in 1:t
        for r in J, s in J
            if r >= s
                continue
            end
            JuMP.@constraints(model, begin
                z[j,r,s] <= x[j,r] + x[j,s]
                z[j,r,s] <= x[j,r] + y[j,r]
                z[j,r,s] <= x[j,s] + y[j,s]
                z[j,r,s] <= y[j,r] + y[j,s]
                z[j,r,s] >= x[j,r] + y[j,s] - 1
                z[j,r,s] >= x[j,s] + y[j,r] - 1
            end)
        end
        for r in J
            JuMP.@constraint(model, x[j,r] + y[j,r] <= 1)
        end
    end

    for r in J, s in J
        if r >= s
            continue
        end
        if E[r,s] == 1
            JuMP.@constraint(model, sum(z[j,r,s] for j in 1:t) == 0)
        else
            JuMP.@constraint(model, sum(z[j,r,s] for j in 1:t) >= 1)
        end
    end

    # Maybe throw in an objctive to maximize spacial logic?
    # JuMP.@objective(model, Min, sum(x) + sum(y))

    stat = JuMP.optimize!(model)

    return termination_status(model) == MOI.OPTIMAL
end

function plot_optimal(E, obstacles, points, lower, upper)
    while upper - lower > 1
        t = (upper + lower) ÷ 2 # integer division
        is_feasible = is_biclique_feasible(E, t)
        
        if is_feasible
            upper = t
        else
            lower = t
        end
    end

    t = upper
    n = size(E)[1]
    J = 1:n

    model = JuMP.Model(Gurobi.Optimizer)
    JuMP.@variable(model, x[1:t,1:n], Bin)
    JuMP.@variable(model, y[1:t,1:n], Bin)
    JuMP.@variable(model, z[1:t,1:n,1:n], Bin)

    for j in 1:t
        for r in J, s in J
            if r >= s
                continue
            end
            JuMP.@constraints(model, begin
                z[j,r,s] <= x[j,r] + x[j,s]
                z[j,r,s] <= x[j,r] + y[j,r]
                z[j,r,s] <= x[j,s] + y[j,s]
                z[j,r,s] <= y[j,r] + y[j,s]
                z[j,r,s] >= x[j,r] + y[j,s] - 1
                z[j,r,s] >= x[j,s] + y[j,r] - 1
            end)
        end
        for r in J
            JuMP.@constraint(model, x[j,r] + y[j,r] <= 1)
        end
    end

    for r in J, s in J
        if r >= s
            continue
        end
        if E[r,s] == 1
            JuMP.@constraint(model, sum(z[j,r,s] for j in 1:t) == 0)
        else
            JuMP.@constraint(model, sum(z[j,r,s] for j in 1:t) >= 1)
        end
    end

    # Maybe throw in an objctive to maximize spacial logic?
    # JuMP.@objective(model, Min, sum(x) + sum(y))

    stat = JuMP.optimize!(model)

    # x[i, j] <- Element i is in left of biclique j
    # y[i, j] <- Element i is in right of biclique j
    if (termination_status(model) == MOI.OPTIMAL)
        for i in 1:size(x)[1]
            Plots.scatter()

            for j in 1:size(x)[2]
                if (value(x[i, j]) == 1)
                    point = points[j]
                    Plots.scatter!([point.first], [point.second], color="red", lims=(-0.1, 1.1), series_annotations=[Plots.text(string(j), :right, 6, "courier")])
                end

                if (value(y[i, j]) == 1)
                    point = points[j]
                    Plots.scatter!([point.first], [point.second], color="blue", lims=(-0.1, 1.1), series_annotations=[Plots.text(string(j), :right, 6, "courier")])
                end
            end


            ClutteredEnvPathOpt.plot_field(obstacles)
            ClutteredEnvPathOpt.plot_lines(obstacles)

            Plots.savefig(string("optimal", i))
        end
    end

    return termination_status(model) == MOI.OPTIMAL
end