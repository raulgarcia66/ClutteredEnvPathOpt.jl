using ClutteredEnvPathOpt
using LinearAlgebra
using Test
using Pipe
using Plots
using JuMP
using Gurobi

################################### solve_deits ########################################
# Create obstacles
num_obs = 3;
seed = 10;
obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)
plot();
ClutteredEnvPathOpt.plot_field(obstacles)
points = ClutteredEnvPathOpt.find_points(obstacles)
ClutteredEnvPathOpt.plot_points(points)
display(plot!())

# Set parameters
N = 30  # number of steps
f1 = [0.0, 0.1, 0.0]  # initial footstep pose 1
f2 = [0.0, 0.0, 0.0]  # initial footstep pose 2
goal = [1.0, 1.0, 0.0]  # goal pose
Q_g = 10 * Matrix{Float64}(I, 3, 3)  # weight between final footstep and goal pose
Q_r = [1 0 0; 0 1 0; 0 0 0] # Matrix{Float64}(I, 3, 3)  # weight between footsteps
q_t = -.05  # weight for trimming unused steps
# TODO: Assure second footstep is within reachability given first

# Optional named arguments
d1 = 0.2
d2 = 0.2
p1 = [0. 0.07]
p2 = [0, -0.27]
delta_x_y_max = 0.1  # max stride norm in space
delta_θ_max = pi/4  # max difference in θ
L = 5  # number of pieces in p.w.l approx. of sine/cosine 

# Compute optimal path
x1, y1, θ1, t1, stats1 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="merged")
x = copy(x1); y = copy(y1); θ = copy(θ1); t = copy(t1);
term_status, r_solve_time, rel_gap, simplex_iters, r_node_count, num_vertices, merged_size, full_size, num_free_faces, num_free_face_ineq, method_used = stats1

x2, y2, θ2, t2, stats2 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="full")
term_status2, r_solve_time2, rel_gap2, simplex_iters2, r_node_count2, num_vertices2, merged_size2, full_size2, num_free_faces2, num_free_face_ineq2, method_used2 = stats2

x3, y3, θ3, t3, stats3 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="bigM")
term_status3, r_solve_time3, rel_gap3, simplex_iters3, r_node_count3, num_vertices3, merged_size3, full_size3, num_free_faces3, num_free_face_ineq3, method_used3 = stats3

seed_range = 1:3
num_obs_range = 1:4
time_diffs = []
for seed in seed_range
    for num_obs in num_obs_range
        println("\nOn test Seed = $seed, Num_Obs = $num_obs")
        obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)
        _,_,_,_,stats = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method = "merged")
        _,_,_,_,stats2 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="full")
        _,_,_,_,stats3 = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="bigM")              
        term_status, solve_time, rel_gap, simplex_iters, node_count, num_vertices, merged_size, full_size, num_free_faces, num_free_face_ineq, method_used = stats
        term_status2, solve_time2, rel_gap2, simplex_iters2, node_count2, num_vertices2, merged_size2, full_size2, num_free_faces2, num_free_face_ineq2, method_used2 = stats2
        term_status3, solve_time3, rel_gap3, simplex_iters3, node_count3, num_vertices3, merged_size3, full_size3, num_free_faces3, num_free_face_ineq3, method_used3 = stats3
        # Merged biclique times vs bigM
        t_diff = solve_time3 - solve_time
        push!(time_diffs, t_diff)
    end
end

# Trim excess steps
num_to_trim = length(filter(tj -> tj > 0.5, t[3:end]))
if num_to_trim % 2 == 0
    x = vcat(x[1:2], x[num_to_trim + 3 : end]);
    y = vcat(y[1:2], y[num_to_trim + 3 : end]);
    θ = vcat(θ[1:2], θ[num_to_trim + 3 : end]);
else
    x = vcat(x[1], x[num_to_trim + 3 : end]);
    y = vcat(y[1], y[num_to_trim + 3 : end]);
    θ = vcat(θ[1], θ[num_to_trim + 3 : end]);
end

num_to_trim2 = length(filter(tj -> tj > 0.5, t2[3:end]))
if num_to_trim2 % 2 == 0
    x2 = vcat(x2[1:2], x2[num_to_trim2 + 3 : end]);
    y2 = vcat(y2[1:2], y2[num_to_trim2 + 3 : end]);
    θ2 = vcat(θ2[1:2], θ2[num_to_trim2 + 3 : end]);
else
    x2 = vcat(x2[1], x2[num_to_trim2 + 3 : end]);
    y2 = vcat(y2[1], y2[num_to_trim2 + 3 : end]);
    θ2 = vcat(θ2[1], θ2[num_to_trim2 + 3 : end]);
end

num_to_trim3 = length(filter(tj -> tj > 0.5, t3[3:end]))
if num_to_trim3 % 2 == 0
    x3 = vcat(x3[1:2], x3[num_to_trim3 + 3 : end]);
    y3 = vcat(y3[1:2], y3[num_to_trim3 + 3 : end]);
    θ3 = vcat(θ3[1:2], θ3[num_to_trim3 + 3 : end]);
else
    x3 = vcat(x3[1], x3[num_to_trim3 + 3 : end]);
    y3 = vcat(y3[1], y3[num_to_trim3 + 3 : end]);
    θ3 = vcat(θ3[1], θ3[num_to_trim3 + 3 : end]);
end

# Plot footstep plan
plot_steps(obstacles, x, y, θ)
plot_steps(obstacles, x2, y2, θ2)
plot_steps(obstacles, x3, y3, θ3)
png("Path Seed $seed Num Obs $num_obs")


# Plot intersections of circles
plot_circles(x, y, θ)
plot_circles(x2, y2, θ2)
plot_circles(x3, y3, θ3)


######################################################################################

num_obs = 3;
seed = 5;
obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs, "Obstacles Seed $seed Num Obs $num_obs", seed = seed)
skeleton = LabeledGraph(g)
all_faces = Set{Vector{Int64}}(union(obstacle_faces, free_faces))

# Plot obstacles
plot(title="Obstacles");
ClutteredEnvPathOpt.plot_field(obstacles)
# points = ClutteredEnvPathOpt.find_points(obstacles)
ClutteredEnvPathOpt.plot_points(points, vertices=skeleton.labels)
# ClutteredEnvPathOpt.plot_lines(obstacles)
# ClutteredEnvPathOpt.plot_borders()
# ClutteredEnvPathOpt.plot_intersections(obstacles);
png("DT Obstacles Seed $seed Num Obs $num_obs")
display(plot!())

# Plot faces
ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, plot_name = "Obstacle Faces", col = "dodgerblue")
png("Obstacle Faces Seed $seed Num Obs $num_obs")

ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "Free Faces")
png("DT Free Faces Seed $seed Num Obs $num_obs")

ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, plot_name = "All Faces", col = "red")
ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "All Faces", new_plot = false)
png("DT All Faces Seed $seed Num Obs $num_obs")

# Plot faces individually
ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, col = "dodgerblue", individually = true)
ClutteredEnvPathOpt.plot_faces(free_faces, points, individually = true)

# Plot edges
ClutteredEnvPathOpt.plot_edges(skeleton, points, plot_name = "Skeleton", vertices=skeleton.labels)
png("DT Skeleton Seed $seed Num Obs $num_obs")

feg_S = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
ClutteredEnvPathOpt.plot_edges(feg_S, points, plot_name = "Finite Element Graph of S",col="black", vertices=feg_S.labels)
png("DT FEG of S Seed $seed Num Obs $num_obs")

feg_obs = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(obstacle_faces))
ClutteredEnvPathOpt.plot_edges(feg_obs, points, plot_name = "Finite Element Graph Obstacle Faces")
png("DT FEG Obstacles Seed $seed Num Obs $num_obs")

feg_all = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces))
ClutteredEnvPathOpt.plot_edges(feg_all, points, plot_name = "Finite Element Graph")
png("DT FEG All Seed $seed Num Obs $num_obs")



# # Duplicate tests
# obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)
# points, mapped, inside_quant = ClutteredEnvPathOpt.find_intersections(obstacles)
# neighbors = ClutteredEnvPathOpt._find_neighbors(points, mapped)
# points = ClutteredEnvPathOpt.find_points(obstacles)

# dup, dup_ind = face_duplicates(obstacle_faces)
# dup, dup_ind = face_duplicates(free_faces)
# dup, dup_ind = face_duplicates(all_faces)
# plot_intersections_indiv(points)
# dup, dup_ind = point_duplicates(points)
# dup, dup_ind = points_in_hs_duplicates(mapped)
# dup, dup_ind = hs_duplicates(mapped)
# points_in_mapped_ordered(obstacles, mapped)
# plot_neighbors(obstacles, neighbors, points)
# list_neighbors(neighbors, points)
# suspect = suspect_neighbors(neighbors)



# cover3 = ClutteredEnvPathOpt._find_biclique_cover_visual(skeleton, free_faces, points)
# tree = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton, free_faces)

# Free Faces Cover
cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces)
(valid_cover, size_diff, missing_edges, missing_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces)), cover)

feg_free = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg_free)
(valid_cover, size_diff, missing_edges, missing_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces)), merged_cover)

ClutteredEnvPathOpt.plot_edges(missing_edges_lg, points, plot_name = "Missing Edges",vertices=missing_edges_lg.labels)
png("Missing Edges Seed $seed Num Obs $num_obs Free")

# file_name = "Biclique Cover Debug Output Seed $seed Num Obs $num_obs Free Unedited.txt"
# f = open(file_name, "w")
# write(f, "Plot 1\n")
# flush(f)
# close(f)
# cover2 = ClutteredEnvPathOpt.find_biclique_cover_debug(skeleton, free_faces, points, obstacles, file_name)

# To produce individual FEGs/Skeletons in the recursion process
(C,A,B), skeleton_ac, faces_ac, skeleton_bc, faces_bc = ClutteredEnvPathOpt.find_biclique_cover_one_iter(skeleton, free_faces)
feg_S_ac = ClutteredEnvPathOpt._find_finite_element_graph(skeleton_ac, ClutteredEnvPathOpt._find_face_pairs(faces_ac))
ClutteredEnvPathOpt.plot_edges(feg_S_ac, points, plot_name = "Finite Element Graph of S_ac",vertices=feg_S_ac.labels,col="black")
# plot!(title="",axis=([], false))
png("DT FEG of S_ac Seed $seed Num Obs $num_obs")
feg_S_bc = ClutteredEnvPathOpt._find_finite_element_graph(skeleton_bc, ClutteredEnvPathOpt._find_face_pairs(faces_bc))
ClutteredEnvPathOpt.plot_edges(feg_S_bc, points, plot_name = "Finite Element Graph of S_bc", vertices=feg_S_bc.labels, col="black")
# plot!(title="",axis=([], false))
png("DT FEG of S_bc Seed $seed Num Obs $num_obs")



# Biclique Cover Validity Tests-----------------------------------------------------------------------
seed_range = 1:100
num_obs_range = 1:4
# fails, tuples = biclique_cover_validity_tests(seed_range, num_obs_range, faces = "all")
fails, tuples = biclique_cover_validity_tests(seed_range, num_obs_range, faces = "free", with_plots=true)

percent_decreases, tuples, c_sizes, mc_sizes = biclique_cover_merger_stats(seed_range, num_obs_range)
percent_decreases_hp, tuples_hp, c_sizes_hp, mc_sizes_hp = biclique_cover_merger_stats(seed_range, num_obs_range, file_name="Biclique Cover Merging Stats HP Partition.txt", partition="hp")
size_diff_cover = c_sizes_hp - c_sizes
size_diff_merged_cover = mc_sizes_hp - mc_sizes

# TODO: WRITE STATS TO FILE

# Biclique validity tester function
function biclique_cover_validity_tests(seed_range, num_obs_range; faces::String="free", with_plots::Bool=false)
    fails = 0
    tuples = []
    for seed = seed_range
        for num_obs = num_obs_range
            println("On test Seed = $seed, Num_Obs = $num_obs")
            obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs,"Obstacles Seed $seed Num Obs $num_obs",seed=seed); # optional seed argument, default is 11
            skeleton = LabeledGraph(g)
            all_faces = union(obstacle_faces, free_faces)
            if faces == "all"
                cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, all_faces);
                feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces));
                merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
                # (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
                (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover);
                println("Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")
            else
                cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces);
                feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces));
                merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
                # (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
                (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover);
                println("Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")
            end
            if size_diff != 0
                fails += 1
                push!(tuples, (seed, num_obs))
                println("Test Failed")

                if with_plots
                    plot(title="Obstacles");
                    ClutteredEnvPathOpt.plot_field(obstacles)
                    # points = ClutteredEnvPathOpt.find_points(obstacles)
                    ClutteredEnvPathOpt.plot_points(points)
                    # ClutteredEnvPathOpt.plot_lines(obstacles)
                    # ClutteredEnvPathOpt.plot_borders()
                    # ClutteredEnvPathOpt.plot_intersections(obstacles);
                    display(plot!())
                    png("Obstacles Seed $seed Num Obs $num_obs")

                    ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "Free Faces")
                    png("Free Faces Seed $seed Num Obs $num_obs")

                    ClutteredEnvPathOpt.plot_edges(skeleton, points, plot_name = "Skeleton")
                    png("Skeleton Seed $seed Num Obs $num_obs")

                    if size_diff > 0
                        ClutteredEnvPathOpt.plot_edges(miss_edges_lg, points, plot_name = "Missing Edges")
                        if faces == "all"
                            png("Missing Edges All Faces Used Seed $seed Num Obs $num_obs")
                        else
                            png("Missing Edges Free Faces Used Seed $seed Num Obs $num_obs")
                        end
                    end
                end
            end
        end
    end
    return (fails, tuples)
end

# Biclique cover merger function for stats
function biclique_cover_merger_stats(seed_range, num_obs_range; file_name="Biclique Cover Merging Stats.txt", partition="CDT")
    # file_name = "Biclique Cover Merging Stats.txt" #Seed Range = $seed_range, Num Obs Range = $num_obs_range.txt"
    f = open(file_name, "w")
    write(f, "Seed Range = $seed_range, Num Obs Range = $num_obs_range\n")
    flush(f)
    close(f)
    
    f = open(file_name, "a")
    write(f, "Seed\tNum_obs\tFull_cover\tMerged_cover\tDecrease\tPercent_decrease\n")
    fails = 0
    tuples = []
    percent_vec = []
    cover_vec = []
    merged_cover_vec = []
    for seed = seed_range
        for num_obs = num_obs_range
            # f = open(file_name, "a")
            # write(f, "\nTest Seed = $seed, Num_Obs = $num_obs\n")
            write(f, "$seed\t$num_obs")
            # flush(f)
            # close(f)

            # println("On test Seed = $seed, Num_Obs = $num_obs")
            obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs,"Obstacles Seed $seed Num Obs $num_obs",seed=seed;partition=partition); # optional seed argument, default is 11
            skeleton = LabeledGraph(g)
            all_faces = union(obstacle_faces, free_faces)
            
            cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces);
            feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces));
            merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
            (valid_cover, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
            (valid_cover_merged, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover);
            # println("Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")

            # if valid_cover && valid_cover_merged
            #     write(f,"Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")
            #     percent = round( ((length(cover)-length(merged_cover)) / length(cover)) * 100, digits = 2)
            #     write(f," | Percent Decrease: $percent\n")
            #     push!(percent_vec, percent)
            # else
            #     fails += 1
            #     push!(tuples, (seed, num_obs))
            #     if !valid_cover_merged
            #         write(f,"Merged cover invalid.\n")
            #     else
            #         write(f,"Original cover invalid.\n")
            #     end
            # end
            if valid_cover && valid_cover_merged
                percent = round( ((length(cover)-length(merged_cover)) / length(cover)) * 100, digits = 2)
                write(f,"\t$(length(cover))\t$(length(merged_cover))\t$(length(cover)-length(merged_cover))\t$percent\n")
                push!(percent_vec, percent)
                push!(cover_vec, length(cover))
                push!(merged_cover_vec, length(merged_cover))
            elseif valid_cover
                percent = round( ((length(cover)-0) / length(cover)) * 100, digits = 2)
                write(f,"\t$(length(cover))\t0\t$(length(cover)-0)\t$percent\n")
                push!(percent_vec, percent)
                push!(cover_vec, -1)
                push!(merged_cover_vec, -1)

                fails += 1
                push!(tuples, (seed, num_obs))
                # if !valid_cover_merged
                #     write(f,"Merged cover invalid.\n")
                # else
                #     write(f,"Original cover invalid.\n")
                # end
            end
        end
    end

    # write(f, "\nTotal Fails: $fails.\nTuples: $tuples")
    flush(f)
    close(f)

    return percent_vec, tuples, cover_vec, merged_cover_vec
end

#------------------------------------------------------------------------------------------------------

# method <- "merged" for the compact biclique cover, "full" for the original biclique cover, "bigM" for big-M constraints
function solve_time_stats(seed_range, num_obs_range; method="merged", file_name="Solve Time Stats.txt")
    # file_name = "Solve Time Stats.txt" #Seed Range = $seed_range, Num Obs Range = $num_obs_range.txt"
    # f = open(file_name, "a")
    # write(f, "Method = $method, Seed Range = $seed_range, Num Obs Range = $num_obs_range\n")
    # flush(f)
    # close(f)

    # Set parameters
    N = 30  # number of steps
    f1 = [0.0, 0.1, 0.0]  # initial footstep pose 1
    f2 = [0.0, 0.0, 0.0]  # initial footstep pose 2
    goal = [1, 1, 0]  # goal pose
    Q_g = 10*Matrix{Float64}(I, 3, 3)  # weight between final footstep and goal pose
    Q_r = [1 0 0; 0 1 0; 0 0 0] # Matrix{Float64}(I, 3, 3)  # weight between footsteps
    q_t = -.05  # weight for trimming unused steps
    # Optional named arguments
    # d1 = 0.2 <- radius of reference foot circle
    # d2 = 0.2 <- radius of moving foot circle
    # p1 = [0, 0.07] <- center of reference foot circle
    # p2 = [0, -0.27] <- center of moving foot circle
    # delta_x_y_max = 0.10  # max stride norm in space
    # delta_θ_max = pi/4  # max difference in θ
    # L = 5  # number of pieces in p.w.l approx. of sine/cosine
    
    f = open(file_name, "a")
    # write(f, "Seed\tNum_obs\tTime\n")
    write(f, "Seed\tNum_obs\tSolve_time\tRel_gap\tSimplex_iterations\tNodes_explored\tNum_vertices\tMerged_size\tFull_size\tNum_free_faces\tNum_free_face_inequalities\n")
    times = zeros(length(seed_range) * length(num_obs_range))
    i = 1
    for seed = seed_range
        for num_obs = num_obs_range
            # f = open(file_name, "a")
            # write(f, "Seed = $seed, Num_Obs = $num_obs\n")
            # flush(f)
            # close(f)
            # println("On test Seed = $seed, Num_Obs = $num_obs")

            # Create obstacles
            obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)

            x, y, θ, t, stats = solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method=method)
            term_status, r_solve_time, rel_gap, simplex_iters, r_node_count, num_vertices, merged_size, full_size, num_free_faces, num_free_face_ineq, method_used = stats

            if method_used == method
                times[i] = r_solve_time # time
            else
                times[i] = -1
            end
            write(f, "$seed\t$num_obs\t$(times[i])\t$rel_gap\t$simplex_iters\t$r_node_count\t$num_vertices\t$merged_size\t$full_size\t$num_free_faces\t$num_free_face_ineq\n")
            flush(f)

            i += 1

            # # Trim excess steps
            # num_to_trim = length(filter(tj -> tj > 0.5, t[3:end]))
            # if num_to_trim % 2 == 0
            #     x = vcat(x[1:2], x[num_to_trim + 3 : end]);
            #     y = vcat(y[1:2], y[num_to_trim + 3 : end]);
            #     θ = vcat(θ[1:2], θ[num_to_trim + 3 : end]);
            # else
            #     x = vcat(x[1], x[num_to_trim + 3 : end]);
            #     y = vcat(y[1], y[num_to_trim + 3 : end]);
            #     θ = vcat(θ[1], θ[num_to_trim + 3 : end]);
            # end
            # plot_steps(obstacles, x, y, θ)
            # plot_circles(x, y, θ, p1=p1, p2=p2)

        end
    end

    # write(f, "\nEND")
    # flush(f)
    close(f)

    return times

end

seed_range = 21:30 # done 1:20
num_obs_range = 1:4
methods = ("merged", "full", "bigM")
times = Dict()
# method = "merged"
for method in methods
    file_name = "Solve Time Stats Method $method.txt"
    f = open(file_name, "a")   # write or append appropriately
    write(f, "Seed Range = $seed_range, Num Obs Range = $num_obs_range\nMethod = $method\n")
    flush(f)
    close(f)
    time = solve_time_stats(seed_range, num_obs_range, method=method, file_name=file_name)
    times[method] = time
end

file_name = "Solve Time Stats Comparison.txt" #Seed Range = $seed_range, Num Obs Range = $num_obs_range.txt"
f = open(file_name, "a")    # write or append appropriately
times_m_b = times["bigM"] .- times["merged"]
times_m_f = times["full"] .- times["merged"]
times_f_b = times["bigM"] .- times["full"]
write(f, "Seed Range = $seed_range, Num Obs Range = $num_obs_range\n")
write(f, "Merged_vs_BigM\tMerged_vs_Full\tFull_vs_BigM\n")
for i = 1:length(times_m_b)
    write(f, "$(times_m_b[i])\t$(times_m_f[i])\t$(times_f_b[i])\n")
end
flush(f)
close(f)


import DelimitedFiles
data_merged = DelimitedFiles.readdlm("Solve Time Stats Method merged.txt", '\t')
matrix_merged = Matrix{Float64}(data_merged[4:end,:])
data_bigM = DelimitedFiles.readdlm("Solve Time Stats Method bigM.txt", '\t')
matrix_bigM = Matrix{Float64}(data_bigM[4:end,:])
header = Vector{String}(data_merged[3,:])

file_name = "Solve Stats Comparison.txt" #Seed Range = $seed_range, Num Obs Range = $num_obs_range.txt"
f = open(file_name, "w")    # write or append appropriately
write(f, "S_T_m\tS_T_b\tS_T_diff\tS_I_m\tS_I_b\tS_I_diff\tN_C_m\tN_C_b\tN_C_diff\n")
S_T_m = matrix_merged[:,3]
S_T_b = matrix_bigM[:,3]
S_T_diff = S_T_b .- S_T_m
S_I_m = matrix_merged[:,5]
S_I_b = matrix_bigM[:,5]
S_I_diff = S_I_b .- S_I_m
N_C_m = matrix_merged[:,6]
N_C_b = matrix_bigM[:,6]
N_C_diff = N_C_b .- N_C_m
for i = 1:length(S_T_m)
    write(f, "$(S_T_m[i])\t$(S_T_b[i])\t$(S_T_diff[i])\t")
    write(f, "$(S_I_m[i])\t$(S_I_b[i])\t$(S_I_diff[i])\t")
    write(f, "$(N_C_m[i])\t$(N_C_b[i])\t$(N_C_diff[i])\n")
end
flush(f)
close(f)

data_summary = DelimitedFiles.readdlm("Solve Stats Comparison.txt", '\t')
matrix_summary = Matrix{Float64}(data_summary[2:end,:])
header_summary = Vector{String}(data_summary[1,:])



# Debugging Functions ---------------------------------------------------------------

# Check if free_faces has any duplicates---------------------------------------------
function face_duplicates(faces)
    dup = 0
    dup_ind = []
    col_faces = collect(faces)
    for i in 1:length(col_faces)
        for j = (i+1):length(col_faces)
            sort_i = sort(col_faces[i])
            sort_j = sort(col_faces[j])
            if sort_i == sort_j
                dup += 1
                push!(dup_ind, (i,j))
            end
        end
    end
    return dup, dup_ind
end

# Check if points has any duplicates---------------------------------------------
function point_duplicates(points)
    dup = 0
    dup_ind = []
    tol = Rational(1.e-13)
    for i in 1:(length(points))
        for j in (i+1):(length(points))
            #if (points[i].first == points[j].first) && (points[i].second == points[j].second)
            if abs(points[i].first - points[j].first) < tol && abs(points[i].second - points[j].second) < tol
                dup += 1
                push!(dup_ind, (i,j))
            end
        end
    end
    # for (i,j) in dup_ind
    #     println("Points $i and $j: \n ($(points[i])) \n ($(points[j]))")
    # end
    return dup, dup_ind
end

# Check if halfspaces in mapped have any duplicated points within them ---------
function points_in_hs_duplicates(mapped)
    dup = 0
    dup_ind = []
    for (l,hs) in enumerate(mapped)
        for i in 1:(length(hs))
            for j in (i+1):(length(hs))
                if (hs[i].first == hs[j].first) && (hs[i].second == hs[j].second)
                    dup += 1
                    push!(dup_ind, (l,i,j))
                end
            end
        end
    end
    return dup, dup_ind
end

# Check if any halfspaces in mapped are duplicates-----------------------------
function hs_duplicates(mapped)
    dup = 0
    dup_ind = []
    for i = 1:length(mapped)
        for j = (i+1):length(mapped)
            if mapped[i] == mapped[j]
                dup += 1
                push!(dup_ind, (i,j))
            end
        end
    end
    return dup, dup_ind
end

# Check to see if points are in order in mapped---------------------------------
function points_in_mapped_ordered(obstacles, mapped)
    for hs in mapped
        plot(legend=false,xlims=(-0.05,1.05),ylims=(-0.05,1.05));
        ClutteredEnvPathOpt.plot_field(obstacles)
        ClutteredEnvPathOpt.plot_lines(obstacles)
        for point in hs
            display(scatter!([point.first], [point.second]))
        end
    end
end

# Plot points individually----------------------------------------------------------
function plot_intersections_indiv(intersections)
    x = map(point -> point[1], intersections)
    y = map(point -> point[2], intersections)
    Plots.plot(legend=false,xlims=(-0.05,1.05),ylims=(-0.05,1.05))
    for i = 1:length(intersections)
        # Plots.plot(legend=false,xlims=(-0.05,1.05),ylims=(-0.05,1.05))
        display(Plots.scatter!([x[i]],[y[i]],title="Point $i:\n$(x[i])\n$(y[i])"))
    end
end

# Check if neighbors are correct----------------------------------------------------
function plot_neighbors(obstacles, neighbors, points)
    for (i,point) in enumerate(points)
        plot();
        ClutteredEnvPathOpt.plot_field(obstacles)
        ClutteredEnvPathOpt.plot_lines(obstacles)
        display(scatter!([point.first],[point.second],title="Point $i"))
        x_i_neighbors = map(node -> points[node].first, neighbors[i])
        y_i_neighbors = map(node -> points[node].second, neighbors[i])
        display(scatter!(x_i_neighbors,y_i_neighbors, title="Point $i. Neighbors $(neighbors[i])"))
        #println("Point $i has neighbors $(neighbors[i]).")
    end
end

function list_neighbors(neighbors, points)
    for i in 1:length(points)
        println("Point $i has neighbors $(neighbors[i]).")
    end
end

function duplicates(n)
    dup = 0
    for i in 1:length(n)
        for j in (i+1):length(n)
            if n[i] == n[j]
                dup += 1
            end
        end
    end
    return dup != 0
end

function suspect_neighbors(neighbors)
    suspect = []
    for n in 1:length(neighbors)
        if length(neighbors[n]) > 4 || duplicates(neighbors[n])
        #if duplicates(neighbors[n])
            push!(suspect, n)
        end
    end
    for point in suspect
        println("Point $point has neighbors $(neighbors[point])")
    end
    return suspect
end

# # Face finder ----------------------------------------------------------------------------
# t_seed = 7
# t_num_obs = 4
# obstacles = ClutteredEnvPathOpt.gen_field(t_num_obs, seed = t_seed) # optional seed
# construct_graph_debug(obstacles)