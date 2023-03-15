using ClutteredEnvPathOpt
using LightGraphs
import Polyhedra
# import GLPK
using Pipe
using Plots
using JuMP, Gurobi, PiecewiseLinearOpt

"""
    get_M_A_b_easy(obstacles)

Computes big-M data by evaluating extreme points of unit cell.
"""
function get_M_A_b_easy(obstacles)
    # TODO: This accepts obstacles and computes points and free faces, where as get_M_A_b accepts points and free_faces
    Ms = []
    As = []
    bs = []
    acc = [0] # accumulator   # try using reduce() instead

    _, points, _, _, free_faces = ClutteredEnvPathOpt.construct_graph(obstacles)

    for face in free_faces
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)   # collect() turns x => y into a vector [x; y]
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

            M = maximum(map(x -> A * x, v.points))

            push!(Ms, M)
            push!(As, A)
            push!(bs, b)
        end
    end
    
    return vcat(Ms...), vcat(As...), vcat(bs...), acc
end

"""
    get_M_A_b(points, free_faces)

Computes big-M data by solving LP over unit cell.
"""
function get_M_A_b(points, free_faces)
    Ms = []
    As = []
    bs = []
    acc = [0] # try using reduce() instead

    # Solve over the unitcell
    u = Polyhedra.convexhull([0//1,0//1],[0//1,1//1],[1//1,0//1],[1//1,1//1])
    unitcell = Polyhedra.polyhedron(u)

    for face in free_faces
        v = Polyhedra.convexhull(map(i -> collect(points[i]), face)...)   # collect() turns x => y into a vector [x; y]
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

            sub_model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
            JuMP.@variable(sub_model, sub_x[1:2])
            JuMP.@constraint(sub_model, sub_x in unitcell)
            JuMP.@objective(sub_model, Max, A * sub_x)
            JuMP.optimize!(sub_model)

            M = JuMP.objective_value(sub_model)

            push!(Ms, M)
            push!(As, A)
            push!(bs, b)
        end
    end

    return vcat(Ms...), vcat(As...), vcat(bs...), acc
end

"""
    plot_circles(x, y, theta; d1=0.20, d2=0.20, p1=[0, 0.07], p2=[0, -0.27])

Plots circles defining reachable region for each footstep onto individual plots.
"""
function plot_circles(x, y, theta; d1=0.20, d2=0.20, p1=[0, 0.07], p2=[0, -0.27])
    angles = LinRange(0,2*pi, 100)
    #plot()
    for i = 2:length(x)
        plot()
        if i % 2 == 1
            circlex1 = d1*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(-p1[1]) - sin(theta[i-1])*(-p1[2]) )
            circley1 = d1*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(-p1[1]) + cos(theta[i-1])*(-p1[2]) )

            circlex2 = d2*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(-p2[1]) - sin(theta[i-1])*(-p2[2]) )
            circley2 = d2*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(-p2[1]) + cos(theta[i-1])*(-p2[2]) )
            #scatter!([x[i-1]; x[i]], [y[i-1]; y[i]])
            scatter!([x[i-1]], [y[i-1]], color = "blue")
            scatter!([x[i]], [y[i]], color = "red")
            quiver!([x[i-1]; x[i]], [y[i-1]; y[i]], 
                quiver=([0.075 * cos(theta[i-1]); 0.075 * cos(theta[i])], 
                [0.075 * sin(theta[i-1]); 0.075 * sin(theta[i])])
                )
        else
            circlex1 = d1*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(p1[1]) - sin(theta[i-1])*(p1[2]) )
            circley1 = d1*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(p1[1]) + cos(theta[i-1])*(p1[2]) )

            circlex2 = d2*cos.(angles) .+ ( x[i-1] + cos(theta[i-1])*(p2[1]) - sin(theta[i-1])*(p2[2]) )
            circley2 = d2*sin.(angles) .+ ( y[i-1] + sin(theta[i-1])*(p2[1]) + cos(theta[i-1])*(p2[2]) )
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
        display(plot!(legend=false, xlims=(-0.05,1.05), ylims=(-0.05,1.05)))   # TODO: Don't hardcode xlims/ylims
    end
    # display(plot!(legend=false, xlims=(-0.05,1.05), ylims=(-0.05,1.05)))
end

"""
    plot_PWL_sin(break_pts)

Plots PWL approximation of sine given breakpoints.
"""
function plot_PWL_sin(break_pts)
    plot(0:0.1:2pi, sin.(0:0.1:2pi), lw=3, title="Piecewise Linear Sine Approximation")

    for i = 2:length(break_pts)
        plot!([break_pts[i-1]; break_pts[i]], [sin(break_pts[i-1]); sin(break_pts[i])], lw=3, color="black")
    end
    display(plot!(legend=false))
end

"""
    plot_PWL_cos(break_pts)

Plots PWL approximation of cosine given breakpoints.
"""
function plot_PWL_cos(break_pts)
    plot(0:0.1:2pi, cos.(0:0.1:2pi), lw=3, title="Piecewise Linear Cosine Approximation")

    for i = 2:length(break_pts)
        plot!([break_pts[i-1]; break_pts[i]], [cos(break_pts[i-1]); cos(break_pts[i])], lw=2, color="black")
    end
    display(plot!(legend=false))
end

# obstacles <- obstacles (list of polyhedra)
# N <- max number of steps (scalar)
# f1 <- initial left foot ([x, y, theta])
# f2 <- initial right foot ([x, y, theta])
# g <- goal pose
# Q_g <- pose weights to goal (3x3 mat)
# Q_r <- pose weights between steps (3x3 mat)
# q_t <- weight on unused steps (scalar)
# method <- "merged" for the compact biclique cover, "full" for the original biclique cover, "bigM" for big-M constraints
# d1 = 0.2 <- radius of reference foot circle
# d2 = 0.2 <- radius of moving foot circle
# p1 = [0, 0.07] <- center of reference foot circle
# p2 = [0, -0.27] <- center of moving foot circle
# delta_x_y_max = 0.1 <- max stride norm in space (no longer used)
# delta_θ_max = pi/4 <- max difference in θ (no longer used)
# relax <- if true, solve as continuous relaxation
"""
    solve_steps(obstacles, N, f1, f2, g, Q_g, Q_r, q_t;
                method="merged", partition="CDT", merge_faces=true, relax=false, 
                MIPGap=0.05, TimeLimit=180, LogFile="", LogToConsole=0,
                d1=0.1, d2=0.1, p1=[0.0, 0.0], p2=[0.0, -0.14])

Compute optimal footstep path for a set of given obstacles and parameters. Adapted from Deits and Tedrake 2014.
Methods:
1) "full" : Independent Branching with full biclique cover computed by our algorithm
1) "merged" : Independent Branching with biclique cover after applying merging procedure to full biclique cover
1) "big-M" : big-M approach
"""
function solve_steps(obstacles, N, f1, f2, g, Q_g, Q_r, q_t; 
    method::String="merged", partition::String="CDT", merge_faces::Bool=false, relax::Bool=false,
    MIPGap=0.05, TimeLimit=180, LogFile="", LogToConsole=0,
    d1=0.1, d2=0.1, p1=[0.0, 0.0], p2=[0.0, -0.14]) #,
    # delta_x_y_max=0.1)

    if partition == "CDT"
        _, points, graph, _, free_faces = ClutteredEnvPathOpt.construct_graph_delaunay(obstacles, merge_faces=merge_faces)
    else
        _, points, graph, _, free_faces = ClutteredEnvPathOpt.construct_graph(obstacles)
    end

    # model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer))
    if relax
        # model = JuMP.Model(Gurobi.Optimizer)
        model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer, "LogFile" => LogFile))
    else
        # TODO: Set MIPGap to large value. Sub optimal solutions still look good
        model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit,"LogFile"=>LogFile)) # "MIPGap"=>MIPGap
        # model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer,"MIPGap"=>MIPGap,"TimeLimit"=>TimeLimit,"LogFile"=>LogFile,"LogToConsole"=>LogToConsole))
        # model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer,"Presolve"=>0, "MIPGap"=>MIPGap,"LogFile"=>LogFile, "TimeLimit"=>TimeLimit)) # "Heuristics"=>0, "Presolve"=>0, "MIPGap"=>.01
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

    # Reachability
    # Breakpoints need to be strategically chosen (fewer is better for problem size)
    s_break_pts = [0, 6pi/16, 10pi/16, 22pi/16, 26pi/16, 2pi]
    c_break_pts = [0, 3pi/16, 14pi/16, 18pi/16, 29pi/16, 2pi]
    # s_break_pts = 0:0.2:2pi  # TODO: Check if maximum strides all the time
    # c_break_pts = 0:0.2:2pi
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

    # Set trimmed steps equal to initial position
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

    # Max step distance. Change to two consecutive (alternating) feet. Distance betweem should be capped at 0.1.
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

    # Footstep location data
    cover = Set{Pair{Set{Int64}, Set{Int64}}}()
    merged_cover = Set{Pair{Set{Int64}, Set{Int64}}}()
    num_free_face_ineq = 0

    if method != "bigM"
        skeleton = LabeledGraph(graph)
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
            JuMP.@constraint(model, x[i] == sum(points[v].first * λ[i, v] for v in 1:J))
            JuMP.@constraint(model, y[i] == sum(points[v].second * λ[i, v] for v in 1:J))
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
            JuMP.@constraint(model, x[i] == sum(points[v].first * λ[i, v] for v in 1:J))
            JuMP.@constraint(model, y[i] == sum(points[v].second * λ[i, v] for v in 1:J))
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
                    JuMP.@constraint(model, A[i, :]' * [x[j]; y[j]] <= b[i] * z[j, r] + M[i] * (1 - z[j, r]))
                end
            end
            JuMP.@constraint(model, sum(z[j, :]) == 1)
        end
    end

    # Solve
    if relax
        relax_integrality(model)
    end
    JuMP.optimize!(model)

    if method == "merged" && valid_cover_merged
        println("\n\nUsed merged cover.")
    elseif method != "bigM" && valid_cover 
        println("\n\nUsed full cover.")
        method = "full"
    else
        println("\n\nUsed big-M constraints.")
        method = "bigM"
    end

    if termination_status(model) == MOI.OPTIMAL
        println("Solution is optimal.\n")        
    elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
        println("Solution is suboptimal due to a time limit, but a primal solution is available.\n")
    else
        error("The model was not solved correctly.\n")
    end

    # stats = (termination_status(model), objective_value(model; result=1), solve_time(model), relative_gap(model), simplex_iterations(model), barrier_iterations(model),
    #         node_count(model), LightGraphs.nv(graph), length(merged_cover), length(cover), length(free_faces), num_free_face_ineq, method)
    stats = (termination_status(model), objective_value(model; result=1), solve_time(model), relative_gap(model),
            node_count(model), LightGraphs.nv(graph), length(merged_cover), length(cover), length(free_faces), num_free_face_ineq, method)

    return value.(x; result=1), value.(y; result=1), value.(θ; result=1), value.(t; result=1), value.(z; result=1), stats
end

"""
    solve_steps(obstacles, N, f1, f2, g, Q_g, Q_r, q_t, x_start, y_start, θ_start, t_start;
                method="merged", partition="CDT", merge_faces=true, relax=false, 
                MIPGap=0.05, TimeLimit=180, LogFile="", LogToConsole=0,
                d1=0.1, d2=0.1, p1=[0.0, 0.0], p2=[0.0, -0.14])

Compute optimal footstep path for a set of given obstacles, parameters, and warm start solutions. Adapted from Deits and Tedrake 2014.
Methods:
1) "full" : Independent Branching with full biclique cover computed by our algorithm
1) "merged" : Independent Branching with biclique cover after applying merging procedure to full biclique cover
1) "big-M" : big-M approach
"""
function solve_steps(obstacles, N, f1, f2, g, Q_g, Q_r, q_t, x_start, y_start, θ_start, t_start;
    method::String="merged", partition::String="CDT", merge_faces::Bool=false, relax::Bool=false,
    MIPGap=0.05, TimeLimit=180, LogFile="", LogToConsole=0,
    d1=0.1, d2=0.1, p1=[0.0, 0.0], p2=[0.0, -0.14]) #,
    # delta_x_y_max=0.1)

    if partition == "CDT"
        _, points, graph, _, free_faces = ClutteredEnvPathOpt.construct_graph_delaunay(obstacles, merge_faces=merge_faces)
    else
        _, points, graph, _, free_faces = ClutteredEnvPathOpt.construct_graph(obstacles)
    end

    # model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer))
    if relax
        # model = JuMP.Model(Gurobi.Optimizer)
        model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer, "LogFile" => LogFile))
    else
        # TODO: Set MIPGap to large value. Sub optimal solutions still look good
        model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit,"LogFile"=>LogFile)) # "MIPGap"=>MIPGap
        # model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer,"MIPGap"=>MIPGap,"TimeLimit"=>TimeLimit,"LogFile"=>LogFile,"LogToConsole"=>LogToConsole))
        # model = JuMP.Model(JuMP.optimizer_with_attributes(Gurobi.Optimizer,"Cuts"=>1,"LogFile"=>LogFile, "TimeLimit"=>TimeLimit)) # "Heuristics"=>0, "Presolve"=>0, "MIPGap"=>.01
    end

    # model has scalar variables x, y, θ, binary variable t
    JuMP.@variable(model, x[1:N])
    JuMP.@variable(model, y[1:N])
    JuMP.@variable(model, θ[1:N])
    JuMP.@variable(model, t[1:N], Bin)

    # Set start values
    set_start_value.(x, x_start)
    set_start_value.(y, y_start)
    set_start_value.(θ, θ_start)
    set_start_value.(t, t_start)

    # Objective
    f = vec([[x[j], y[j], θ[j]] for j in 1:N])
    # f_no_theta = vec([[x[j], y[j]] for j in 1:N])
    JuMP.@objective(
        model,
        Min,
        ((f[N] - g)' * Q_g * (f[N] - g)) + sum(q_t * t) + sum((f[j + 1] - f[j])' * Q_r * (f[j + 1] - f[j]) for j in 1:(N-1))
    )

    # Reachability
    # Breakpoints need to be strategically chosen (fewer is better for problem size)
    s_break_pts = [0, 6pi/16, 10pi/16, 22pi/16, 26pi/16, 2pi]
    c_break_pts = [0, 3pi/16, 14pi/16, 18pi/16, 29pi/16, 2pi]
    # s_break_pts = 0:0.2:2pi  # TODO: Check if maximum strides all the time
    # c_break_pts = 0:0.2:2pi
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

    # Set trimmed steps equal to initial position
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

    # Max step distance. Change to two consecutive (alternating) feet. Distance betweem should be capped at 0.1.
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

    # Footstep location data
    cover = Set{Pair{Set{Int64}, Set{Int64}}}()
    merged_cover = Set{Pair{Set{Int64}, Set{Int64}}}()
    num_free_face_ineq = 0

    if method != "bigM"
        skeleton = LabeledGraph(graph)
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
            JuMP.@constraint(model, x[i] == sum(points[v].first * λ[i, v] for v in 1:J))
            JuMP.@constraint(model, y[i] == sum(points[v].second * λ[i, v] for v in 1:J))
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
            JuMP.@constraint(model, x[i] == sum(points[v].first * λ[i, v] for v in 1:J))
            JuMP.@constraint(model, y[i] == sum(points[v].second * λ[i, v] for v in 1:J))
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
                    JuMP.@constraint(model, A[i, :]' * [x[j]; y[j]] <= b[i] * z[j, r] + M[i] * (1 - z[j, r]))
                end
            end
            JuMP.@constraint(model, sum(z[j, :]) == 1)
        end
    end

    # Solve
    if relax
        relax_integrality(model)
    end
    JuMP.optimize!(model)

    if method == "merged" && valid_cover_merged
        println("\n\nUsed merged cover.")
    elseif method != "bigM" && valid_cover 
        println("\n\nUsed full cover.")
        method = "full"
    else
        println("\n\nUsed big-M constraints.")
        method = "bigM"
    end

    if termination_status(model) == MOI.OPTIMAL
        println("Solution is optimal.\n")        
    elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
        println("Solution is suboptimal due to a time limit, but a primal solution is available.\n")
    else
        error("The model was not solved correctly.\n")
    end

    # stats = (termination_status(model), objective_value(model; result=1), solve_time(model), relative_gap(model), simplex_iterations(model), barrier_iterations(model),
    #         node_count(model), LightGraphs.nv(graph), length(merged_cover), length(cover), length(free_faces), num_free_face_ineq, method)
    stats = (termination_status(model), objective_value(model; result=1), solve_time(model), relative_gap(model),
            node_count(model), LightGraphs.nv(graph), length(merged_cover), length(cover), length(free_faces), num_free_face_ineq, method)

    return value.(x; result=1), value.(y; result=1), value.(θ; result=1), value.(t; result=1), value.(z; result=1), stats
end