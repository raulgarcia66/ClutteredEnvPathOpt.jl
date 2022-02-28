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
seed = 5;
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
x1, y1, θ1, t1, stats1 = ClutteredEnvPathOpt.solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="merged")
x = copy(x1); y = copy(y1); θ = copy(θ1); t = copy(t1);
term_status, r_solve_time, rel_gap, simplex_iters, r_node_count, num_vertices, merged_size, full_size, num_free_faces, num_free_face_ineq, method_used = stats1

x2, y2, θ2, t2, stats2 = ClutteredEnvPathOpt.solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="full")
term_status2, r_solve_time2, rel_gap2, simplex_iters2, r_node_count2, num_vertices2, merged_size2, full_size2, num_free_faces2, num_free_face_ineq2, method_used2 = stats2

x3, y3, θ3, t3, stats3 = ClutteredEnvPathOpt.solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="bigM")
term_status3, r_solve_time3, rel_gap3, simplex_iters3, r_node_count3, num_vertices3, merged_size3, full_size3, num_free_faces3, num_free_face_ineq3, method_used3 = stats3

seed_range = 1:3
num_obs_range = 1:4
time_diffs = []
for seed in seed_range
    for num_obs in num_obs_range
        println("\nOn test Seed = $seed, Num_Obs = $num_obs")
        obstacles = ClutteredEnvPathOpt.gen_field(num_obs, seed = seed)
        _,_,_,_,stats = ClutteredEnvPathOpt.solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method = "merged")
        _,_,_,_,stats2 = ClutteredEnvPathOpt.solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="full")
        _,_,_,_,stats3 = ClutteredEnvPathOpt.solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method="bigM")              
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
ClutteredEnvPathOpt.plot_steps(obstacles, x, y, θ)
ClutteredEnvPathOpt.plot_steps(obstacles, x2, y2, θ2)
ClutteredEnvPathOpt.plot_steps(obstacles, x3, y3, θ3)
png("Path Seed $seed Num Obs $num_obs")


# Plot intersections of circles
ClutteredEnvPathOpt.plot_circles(x, y, θ)
ClutteredEnvPathOpt.plot_circles(x2, y2, θ2)
ClutteredEnvPathOpt.plot_circles(x3, y3, θ3)


######################################################################################

num_obs = 4;
seed = 3;
merge_faces = false;
partition = "CDT";
obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs, "Obstacles Seed $seed Num Obs $num_obs", seed=seed, partition=partition, merge_faces=merge_faces)
skeleton = LabeledGraph(g)
all_faces = Set{Vector{Int64}}(union(obstacle_faces, free_faces))

# Plot obstacles
plot(title="Obstacles");
ClutteredEnvPathOpt.plot_field(obstacles)
# points = ClutteredEnvPathOpt.find_points(obstacles)
# ClutteredEnvPathOpt.plot_points(points, vertices=skeleton.labels)
# ClutteredEnvPathOpt.plot_lines(obstacles)
# ClutteredEnvPathOpt.plot_borders()
# ClutteredEnvPathOpt.plot_intersections(obstacles);
# plot!(title="",axis=([], false))
png("$partition Obstacles Seed $seed Num Obs $num_obs")
# png("$partition Obstacles Seed $seed Num Obs $num_obs No Axis")
ClutteredEnvPathOpt.plot_points(points, vertices=skeleton.labels)
png("$partition Obstacles with Points Seed $seed Num Obs $num_obs")
# png("$partition Obstacles with Points Seed $seed Num Obs $num_obs No Axis")
display(plot!())

# Plot faces
ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, plot_name = "Obstacle Faces", col = "dodgerblue")
# plot!(title="",axis=([], false))
png("Obstacle Faces Seed $seed Num Obs $num_obs")

ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "Free Faces")
# plot!(title="",axis=([], false))
png("$partition Free Faces Seed $seed Num Obs $num_obs Merged Faces $merge_faces")
# png("$partition Free Faces Seed $seed Num Obs $num_obs Merged Faces $merge_faces No Axis")

ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, plot_name = "All Faces", col = "red")
ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "All Faces", new_plot = false)
# plot!(title="",axis=([], false))
png("$partition All Faces Seed $seed Num Obs $num_obs")

# Plot faces individually
ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, col = "dodgerblue", individually = true)
ClutteredEnvPathOpt.plot_faces(free_faces, points, individually = true)

# Plot edges
ClutteredEnvPathOpt.plot_edges(skeleton, points, plot_name = "Skeleton", col="black", vertices=skeleton.labels)
# plot!(title="",axis=([], false))
png("$partition Skeleton Seed $seed Num Obs $num_obs Merged Faces $merge_faces")
# png("$partition Skeleton Seed $seed Num Obs $num_obs Merged Faces $merge_faces No Axis")

feg_S = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
ClutteredEnvPathOpt.plot_edges(feg_S, points, plot_name = "Finite Element Graph of S",col="black", vertices=feg_S.labels)
# plot!(title="",axis=([], false))
png("$partition FEG of S Seed $seed Num Obs $num_obs Merged Faces $merge_faces")
# png("$partition FEG of S Seed $seed Num Obs $num_obs Merged Faces $merge_faces No Axis")

# feg_obs = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(obstacle_faces))
# ClutteredEnvPathOpt.plot_edges(feg_obs, points, plot_name = "Finite Element Graph Obstacle Faces")
# # plot!(title="",axis=([], false))
# png("$partition FEG Obstacles Seed $seed Num Obs $num_obs")

# feg_all = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces))
# ClutteredEnvPathOpt.plot_edges(feg_all, points, plot_name = "Finite Element Graph")
# # plot!(title="",axis=([], false))
# png("$partition FEG All Seed $seed Num Obs $num_obs")



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
feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
(valid_cover, size_diff, missing_edges, missing_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover)

merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
(valid_cover, size_diff, missing_edges, missing_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover)


ClutteredEnvPathOpt.plot_edges(missing_edges_lg, points, plot_name = "Missing Edges",vertices=missing_edges_lg.labels)
png("$partitition Missing Edges Seed $seed Num Obs $num_obs Merged Faces $merge_faces")

# file_name = "Biclique Cover Debug Output Seed $seed Num Obs $num_obs Free Unedited.txt"
# f = open(file_name, "w")
# write(f, "Plot 1\n")
# flush(f)
# close(f)
# cover2 = ClutteredEnvPathOpt.find_biclique_cover_debug(skeleton, free_faces, points, obstacles, file_name)

# To produce individual FEGs/Skeletons in the recursion process
(C,A,B), skeleton_ac, faces_ac, skeleton_bc, faces_bc = ClutteredEnvPathOpt.find_biclique_cover_one_iter(skeleton, free_faces)
feg_S_ac = ClutteredEnvPathOpt._find_finite_element_graph(skeleton_ac, ClutteredEnvPathOpt._find_face_pairs(faces_ac))
feg_S_bc = ClutteredEnvPathOpt._find_finite_element_graph(skeleton_bc, ClutteredEnvPathOpt._find_face_pairs(faces_bc))

ClutteredEnvPathOpt.plot_edges(skeleton_ac, points, plot_name = "Skeleton S_ac",vertices=skeleton_ac.labels,col="black")
# plot!(title="",axis=([], false))
png("$partition Skeleton of S_ac Seed $seed Num Obs $num_obs Merged Faces $merge_faces")
# png("$partition Skeleton of S_ac Seed $seed Num Obs $num_obs Merged Faces $merge_faces No Axis")
ClutteredEnvPathOpt.plot_edges(feg_S_ac, points, plot_name = "Finite Element Graph of S_ac",vertices=feg_S_ac.labels,col="black")
# plot!(title="",axis=([], false))
png("$partition FEG of S_ac Seed $seed Num Obs $num_obs Merged Faces $merge_faces")
# png("$partition FEG of S_ac Seed $seed Num Obs $num_obs Merged Faces $merge_faces No Axis")

ClutteredEnvPathOpt.plot_edges(skeleton_bc, points, plot_name = "Skeleton S_bc",vertices=skeleton_bc.labels,col="black")
# plot!(title="",axis=([], false))
png("$partition Skeleton of S_bc Seed $seed Num Obs $num_obs Merged Faces $merge_faces")
# png("$partition Skeleton of S_bc Seed $seed Num Obs $num_obs Merged Faces $merge_faces No Axis")
ClutteredEnvPathOpt.plot_edges(feg_S_bc, points, plot_name = "Finite Element Graph of S_bc",vertices=feg_S_bc.labels,col="black")
# plot!(title="",axis=([], false))
png("$partition FEG of S_bc Seed $seed Num Obs $num_obs Merged Faces $merge_faces")
# png("$partition FEG of S_bc Seed $seed Num Obs $num_obs Merged Faces $merge_faces No Axis")



########################### Biclique Cover Validity Tests #############################
seed_range = 1:100
num_obs_range = 1:4
partition = "CDT"
merge_faces = true
# fail_tuples = biclique_cover_validity_tests(seed_range, num_obs_range, faces = "all", partition=partition)
fail_tuples = biclique_cover_validity_tests(seed_range, num_obs_range, faces = "free", with_plots=true, partition=partition)

# Biclique validity tester function
function biclique_cover_validity_tests(seed_range, num_obs_range; faces::String="free", with_plots::Bool=false,partition="CDT",merge_faces::Bool=true)
    # fails = 0
    tuples = []
    for seed = seed_range
        for num_obs = num_obs_range
            println("On test Seed = $seed, Num_Obs = $num_obs")
            obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs,"Obstacles Seed $seed Num Obs $num_obs",seed=seed, partition=partition, merge_faces=merge_faces)
            skeleton = LabeledGraph(g)
            all_faces = Set{Vector{Int64}}(union(obstacle_faces, free_faces))
            if faces == "all"
                cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, all_faces);
                feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces));
                merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
                # (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
                (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover);
                # println("Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")
            else
                cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces);
                feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces));
                merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
                # (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
                (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover);
                # println("Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")
            end
            if size_diff != 0
                # fails += 1
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
                    png("$partition Obstacles Seed $seed Num Obs $num_obs")

                    ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "Free Faces")
                    png("$partition Free Faces Seed $seed Num Obs $num_obs")

                    ClutteredEnvPathOpt.plot_edges(skeleton, points, plot_name = "Skeleton")
                    png("$partition Skeleton Seed $seed Num Obs $num_obs")

                    if size_diff > 0
                        ClutteredEnvPathOpt.plot_edges(miss_edges_lg, points, plot_name = "Missing Edges")
                        if faces == "all"
                            png("$partition Missing Edges All Faces Used Seed $seed Num Obs $num_obs")
                        else
                            png("$partition Missing Edges Free Faces Used Seed $seed Num Obs $num_obs")
                        end
                    end
                end
            end
        end
    end
    return tuples
end

######################### Biclique Cover Merging Stats #######################

file_name = "Biclique Cover Merging Stats Partition $partition.txt"
file_name_dt = "Biclique Cover Merging Stats Partition CDT.txt"
file_name_hp = "Biclique Cover Merging Stats Partition HP.txt"
percent_decreases_dt, tuples_dt, c_sizes_dt, mc_sizes_dt, num_free_faces_dt = biclique_cover_merger_stats(seed_range, num_obs_range, file_name=file_name_dt, partition="CDT")
percent_decreases_hp, tuples_hp, c_sizes_hp, mc_sizes_hp, num_free_faces_hp = biclique_cover_merger_stats(seed_range, num_obs_range, file_name=file_name_hp, partition="HP")
size_diff_cover = c_sizes_hp - c_sizes_dt
size_diff_merged_cover = mc_sizes_hp - mc_sizes_dt
size_diff_free_faces = num_free_faces_hp - num_free_faces_dt

maximum(size_diff_free_faces) # 68
minimum(size_diff_free_faces) # 0
Statistics.mean(size_diff_free_faces) # 17.79
count(t-> t > 0, size_diff_free_faces) # 373
count(t-> t < 0, size_diff_free_faces) # 0
count(t-> t == 0, size_diff_free_faces) # 27

maximum(size_diff_merged_cover) # 16
minimum(size_diff_merged_cover) # -1
Statistics.mean(size_diff_merged_cover) # 4.975
count(t-> t > 0, size_diff_merged_cover) # 353
count(t-> t < 0, size_diff_merged_cover) # 9
count(t-> t == 0, size_diff_merged_cover) # 38

# Biclique cover merger function for stats
function biclique_cover_merger_stats(seed_range, num_obs_range; file_name="Biclique Cover Merging Stats.txt", partition="CDT", merge_faces=true)
    # file_name = "Biclique Cover Merging Stats.txt" #Seed Range = $seed_range, Num Obs Range = $num_obs_range.txt"
    f = open(file_name, "w")
    write(f, "Seed Range = $seed_range, Num Obs Range = $num_obs_range\n")
    flush(f)
    close(f)
    
    f = open(file_name, "a")
    write(f, "Seed\tNum_obs\tNum_free_faces\tFull_cover\tMerged_cover\tDecrease\tPercent_decrease\n")
    tuples = []
    num_free_faces_vec = []
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
            obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs,"Obstacles Seed $seed Num Obs $num_obs",seed=seed,partition=partition,merge_faces=merge_faces);
            skeleton = LabeledGraph(g)
            # all_faces = Set{Vector{Int64}}(union(obstacle_faces, free_faces))
            write(f, "\t$(length(free_faces))")
            push!(num_free_faces_vec, (length(free_faces)))
            
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

                push!(tuples, (seed, num_obs))
                # if !valid_cover_merged
                #     write(f,"Merged cover invalid.\n")
                # else
                #     write(f,"Original cover invalid.\n")
                # end
            end
        end
    end

    # write(f, "\nTotal Fails: $length(tuples).\nTuples: $tuples")
    flush(f)
    close(f)

    return percent_vec, tuples, cover_vec, merged_cover_vec, num_free_faces_vec
end

# TODO: Add number of points, ratio (percentage) of num_free_faces for HP vs DT
# file_name = "Biclique Cover Merging and Free Faces Comparison.txt"
# f = open(file_name, "w")
# write(f, "Title\n")
# write(f, "Header1\tHeader2\n")
for i = 1:length(percent_decreases_dt)
    # percent decreases etc
    write()
end
# flush(f)
# close(f)

# TODO: WRITE STATS TO FILE
# function size_stats()

# end

######################### Solve time function ###########################

# method <- "merged" for the compact biclique cover, "full" for the original biclique cover, "bigM" for big-M constraints
function solve_time_stats(seed_range, num_obs_range; file_name="Solve Time Stats.txt", method="merged", partition="CDT", merge_faces=true)
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
    write(f, "Seed\tNum_obs\tSolve_time\tRel_gap\tSimplex_iterations\tNodes_explored\tNum_vertices\tBC_merged_size\tBC_full_size\tNum_free_faces\tNum_free_face_inequalities\n")
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

            x, y, θ, t, stats = ClutteredEnvPathOpt.solve_deits(obstacles, N, f1, f2, goal, Q_g, Q_r, q_t, method=method, partition=partition, merge_faces=merge_faces)
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

########################### Solve Time Stats ###############################

seed_start = 26 #51
seed_end = 50  #100
seed_range = seed_start:seed_end
num_obs = 4
num_obs_range = num_obs:num_obs
# partitions = ["CDT", "HP"]
partitions = ["CDT"]
# merge_faces = [true, false]
merge_faces = [true]
# merge_faces = [false]
# methods = ["merged", "full", "bigM"]
methods = ["merged"]
# methods = ["full"]
# methods = ["bigM"]
# times = Dict()
for partition in partitions
    if partition == "CDT"
        for merge_face in merge_faces
            for method in methods
                file_name = "Solve Time Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method $method Partition $partition Merge Face $merge_face.txt"
                f = open(file_name, "w")   # write or append appropriately
                write(f, "Seed Range = $seed_range, Num Obs Range = $num_obs_range\nMethod = $method\tPartition = $partition\tMerge Face = $merge_face\n")
                flush(f)
                close(f)
                _ = solve_time_stats(seed_range, num_obs_range, file_name=file_name, method=method, partition=partition, merge_faces=merge_face)
                # times[method] = time
            end
        end
    else
        for method in methods
            file_name = "Solve Time Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method $method Partition $partition.txt"
            f = open(file_name, "w")   # write or append appropriately
            write(f, "Seed Range = $seed_range, Num Obs Range = $num_obs_range\nMethod = $method\tPartition = $partition\n")
            flush(f)
            close(f)
            _ = solve_time_stats(seed_range, num_obs_range, file_name=file_name, method=method, partition=partition)
            # times[method] = time
        end
    end
end

# Import CDT data into matrices
import DelimitedFiles

seed_start = 1
seed_end = 50
num_obs = 1

# CDT, Merged Faces, Merged BC
filename = "Solve Time Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method merged Partition CDT Merge Face true.txt"
CDT_MF_MBC_all_data = DelimitedFiles.readdlm(filename, '\t')
# CDT, Merged Faces, Full BC
filename = "Solve Time Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method full Partition CDT Merge Face true.txt"
CDT_MF_FBC_all_data = DelimitedFiles.readdlm(filename, '\t')
# CDT, Merged Faces, BigM
filename = "Solve Time Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method bigM Partition CDT Merge Face true.txt"
CDT_MF_BigM_all_data = DelimitedFiles.readdlm(filename, '\t')
# CDT, All Faces, Merged BC
filename = "Solve Time Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method merged Partition CDT Merge Face false.txt"
CDT_AF_MBC_all_data = DelimitedFiles.readdlm(filename, '\t')
# CDT, All Faces, Full BC
filename = "Solve Time Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method full Partition CDT Merge Face false.txt"
CDT_AF_FBC_all_data = DelimitedFiles.readdlm(filename, '\t')
# CDT, All Faces, BigM
filename = "Solve Time Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method bigM Partition CDT Merge Face false.txt"
CDT_AF_BigM_all_data = DelimitedFiles.readdlm(filename, '\t')

header = Vector{String}(CDT_MF_MBC_all_data[3,:])
# column indices: 1,2,3 (10,11,12) = seed, num_obs, solve_time (MBC_size, FBC_size, Num_Free_Faces)
CDT_MF_MBC = CDT_MF_MBC_all_data[4:end,3]
CDT_MF_FBC = CDT_MF_FBC_all_data[4:end,3]
CDT_MF_BigM = CDT_MF_BigM_all_data[4:end,3]
CDT_AF_MBC = CDT_AF_MBC_all_data[4:end,3]
CDT_AF_FBC = CDT_AF_FBC_all_data[4:end,3]
CDT_AF_BigM = CDT_AF_BigM_all_data[4:end,3]
test_names = ["MF_MBC","MF_FBC","MF_BigM","AF_MBC","AF_FBC","AF_BigM"]
test_matrix = hcat(CDT_MF_MBC, CDT_MF_FBC, CDT_MF_BigM, CDT_AF_MBC, CDT_AF_FBC, CDT_AF_BigM)

# x1 = CDT_MF_MBC_all_data[4:end,8:10]
# x2 = CDT_MF_FBC_all_data[4:end,8:10]
# x3 = CDT_MF_BigM_all_data[4:end,8:10]
# x4 = CDT_AF_MBC_all_data[4:end,8:10]
# x5 = CDT_AF_FBC_all_data[4:end,8:10]
# x6 = CDT_AF_BigM_all_data[4:end,8:10]

winners = Vector{String}(undef,50)
for i=1:50
    index = argmin(test_matrix[i,:])
    winners[i] = test_names[index]
end
final_matrix = hcat(test_matrix, winners)

winner_count = zeros(Int,6)
for (i,test) in enumerate(test_names)
    winner_count[i] = count(w -> w == test, winners)
end
# For seed 1:50, num obs 1: [10, 10, 19, 3, 1, 7], meaning bigM fastest, full and merged BC were same, and merged faces usually faster
# For seed 1:50, num obs 2: [10, 7, 24, 4, 0, 5], meaning bigM fastest, then merged BC, fullBC, and merged faces usually faster


# Write summary of important data
file_name = "Solve Time Stats Comparison CDT Seed Range $seed_start to $seed_end Num Obs $num_obs.txt"
f = open(file_name, "w")    # write or append appropriately
write(f, "Partition CDT, Seed Range = $seed_start:$seed_end, Num Obs = $num_obs\n")
write(f,"Seed\tNum_Obs\t")
for i = 1:(length(test_names))
    write(f, "$(test_names[i])\t")
end
write(f, "Fastest Method\n")
flush(f)

for i = 1:size(test_matrix, 1)
    write(f,"$i\t1\t")   # seed and num_obs
    for j = 1:size(test_matrix, 2)
        write(f, "$(test_matrix[i,j])\t")
    end
    write(f, "$(winners[i])\n")
end
flush(f)
close(f)



# To load new data summary into matrix
file_name = "Solve Time Stats Comparison CDT Seed Range $seed_start to $seed_end Num Obs $num_obs.txt"
data_summary = DelimitedFiles.readdlm(file_name, '\t')
matrix_summary = data_summary[3:end,:]
num_matrix_summary = Matrix{Float64}(matrix_summary[:,3:end-1])
header_summary = Vector{String}(data_summary[2,:])
winners = data_summary[3:end,end]

winner_count = zeros(Int,6)
for (i,test) in enumerate(test_names)
    winner_count[i] = count(w -> w == test, winners)
end

# mergedBC vs bigM (1 vs 2) (Merged Faces both)
MBC_wins = 0
FBC_wins = 0
BigM_wins = 0
for i = 1:50
    # Don't think this is exactly correct for counting
    if num_matrix_summary[i,1] < num_matrix_summary[i,3]
        MBC_wins += 1     
    elseif num_matrix_summary[i,2] < num_matrix_summary[i,3]
        FBC_wins += 1  
    else
        BigM_wins +=1 
    end
end
MBC_wins
FBC_wins
BigM_wins
# Num Obs 1: 24 vs 5 vs 21
# Num Obs 2:  14 vs 7 vs 29


# # Previous data ------------------------------------------------------------------
# import DelimitedFiles
# data_merged = DelimitedFiles.readdlm("Solve Time Stats Method merged.txt", '\t')
# matrix_merged = Matrix{Float64}(data_merged[4:end,:])
# data_bigM = DelimitedFiles.readdlm("Solve Time Stats Method bigM.txt", '\t')
# matrix_bigM = Matrix{Float64}(data_bigM[4:end,:])
# header = Vector{String}(data_merged[3,:])

# file_name = "Solve Stats Comparison.txt" #Seed Range = $seed_range, Num Obs Range = $num_obs_range.txt"
# f = open(file_name, "w")    # write or append appropriately
# write(f, "S_T_m\tS_T_b\tS_T_diff\tS_I_m\tS_I_b\tS_I_diff\tN_C_m\tN_C_b\tN_C_diff\n")
# S_T_m = matrix_merged[:,3]
# S_T_b = matrix_bigM[:,3]
# S_T_diff = S_T_b .- S_T_m
# S_I_m = matrix_merged[:,5]
# S_I_b = matrix_bigM[:,5]
# S_I_diff = S_I_b .- S_I_m
# N_C_m = matrix_merged[:,6]
# N_C_b = matrix_bigM[:,6]
# N_C_diff = N_C_b .- N_C_m
# for i = 1:length(S_T_m)
#     write(f, "$(S_T_m[i])\t$(S_T_b[i])\t$(S_T_diff[i])\t")
#     write(f, "$(S_I_m[i])\t$(S_I_b[i])\t$(S_I_diff[i])\t")
#     write(f, "$(N_C_m[i])\t$(N_C_b[i])\t$(N_C_diff[i])\n")
# end
# flush(f)
# close(f)

# data_summary = DelimitedFiles.readdlm("Solve Stats Comparison.txt", '\t')
# matrix_summary = Matrix{Float64}(data_summary[2:end,:])
# header_summary = Vector{String}(data_summary[1,:])



#######################################################################################
################################ Debugging Functions ##################################

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
