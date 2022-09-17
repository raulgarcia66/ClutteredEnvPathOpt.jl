using ClutteredEnvPathOpt
using LinearAlgebra
using Test
using Pipe
using Plots
using JuMP
# using Gurobi

# import GLPK
using LightGraphs
import Polyhedra
# using PiecewiseLinearOpt

######################################################################################
########################### Save Obstacle Data from Random ###########################
## THIS IS DONE

# one_obs_seeds = Set([24,30,33,42])
# two_obs_seeds = Set([4,6,8,12,14,15,19,20,21,22,25,28,34,36,37,39,41,46,47,49,51,53,63,64,65,72,73,75,76,77,79,81,82,83,84,85,87,88,90,91,94,98,99])
# three_obs_seeds = setdiff( setdiff( Set(1:100), two_obs_seeds), one_obs_seeds)

# # Original seeds before renaming some meh to good (updated in "Obstacle Selection.jl")
# good_seeds = Set([100,98,96,94,93,92,90,89,86,83,80,78,77,75,73,70,69,68,67,66,65,64,63,62,58,57,56,55,54,53,51,50,48,47,46,45,43,42,39,37,36,35,34,33,32,30,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,7,6,5,4,2,1])
# meh_seeds = Set([99,97,95,91,88,87,84,82,76,72,71,60,59,52,49,44,41,38,28,27,8,3])
# bad_seeds = Set([85,81,79,74,61,40,31,29])

# good_three_obs_seeds = intersect(three_obs_seeds, good_seeds)
# good_two_obs_seeds = intersect(two_obs_seeds, good_seeds)
# good_one_obs_seeds = intersect(one_obs_seeds, good_seeds)

# num_obs = 3
# merge_faces = false;
# partition = "CDT";
# # seeds = setdiff(good_two_obs_seeds, Set([33, 21, 63, 65, 90, 19,4,46,34]))
# for seed in seeds
#     # seed = 5;
#     obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs, "Obstacles Seed $seed Num Obs $num_obs", seed=seed, partition=partition, merge_faces=merge_faces)
#     # Line below is for when you already have a collection of obstacles
#     # obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(obstacles, "Obstacles Seed $seed Num Obs $num_obs", partition=partition, merge_faces=merge_faces)
#     # skeleton = LabeledGraph(g)
#     # all_faces = Set{Vector{Int64}}(union(obstacle_faces, free_faces))
#     plot(title="Obstacles Seed $seed");
#     ClutteredEnvPathOpt.plot_field(obstacles)
#     # png("./test/obstacle files/Seed $seed")
#     x = [i//21 for i = 0:21]
#     for j = 0:21
#         y = [j//21 for i = 0:21]
#         scatter!(x,y)
#     end
#     display(plot!())
#     # TODO: Save a plot of the grid and some plots of obstacles on top of the grid

#     # # Code for saving points to files
#     # file_name = "./test/obstacle files/Seed $seed.txt"
#     # f = open(file_name, "w")   # write or append appropriately
#     # for ob in obstacles
#     #     pts = Polyhedra.points(ob)
#     #     for pt in pts
#     #         # println("$(pt[1]), $(pt[2])")
#     #         write(f, "$(pt[1]), $(pt[2])\n")
#     #     end
#     #     write(f, "END\n")
#     #     flush(f)
#     #     # close(f)
#     # end
#     # # flush(f)
#     # close(f)
# end


# # Manually save obstacle data to files
# seed = 103
# v = Polyhedra.convexhull([4//7,2//7], [3//7, 1//3], [3//7, 11//21], [10//21, 2//3], [4//7, 5//7])
# poly = Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
# plot(title="Obstacles Seed $seed")
# ClutteredEnvPathOpt.plot_field([obstacles[1], obstacles[2], poly])
# display(plot!())
# png("./test/obstacle files/Seed $seed")

# file_name = "./test/obstacle files/Seed $seed.txt"
# f = open(file_name, "w")   # write or append appropriately
# for ob in obstacles
#     pts = Polyhedra.points(ob)
#     for pt in pts
#         # println("$(pt[1]), $(pt[2])")
#         write(f, "$(pt[1]), $(pt[2])\n")
#     end
#     write(f, "END\n")
#     flush(f)
#     # close(f)
# end
# # flush(f)
# close(f)

######################################################################################
################################ Inspect Instances ###################################

## Create from Random package
# num_obs = 3;
# seed = 71;
# merge_faces = false;
# partition = "CDT";
# obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs, "Obstacles Seed $seed Num Obs $num_obs", seed=seed, partition=partition, merge_faces=merge_faces)

## Load from files

function obs_from_file(seed, num_obs, file_name; display_plot=true, save_plot=false)
    obstacles = []
    obs_counter = 0

    points = Vector{Rational{Int}}[]
    for line in eachline(file_name)
        if line != "END"
            pt = split(line, ",")
            # println("$(pt[1]), $(pt[2])")
            pts_rational = map(coord -> parse.(Int, split(coord, "//")) , pt)  # add '/' as well?
            # println("$(typeof(pts_rational))")
            # println("$(Rational(pts_rational[1][1], pts_rational[1][2]))")
            # println("$(typeof([Rational(pts_rational[1][1], pts_rational[1][2]) ; Rational(pts_rational[2][1], pts_rational[2][2])]))")
            push!(points, [Rational(pts_rational[1][1], pts_rational[1][2]) ; Rational(pts_rational[2][1], pts_rational[2][2])])
        else
            # v = Polyhedra.convexhull(points...)
            # poly = Polyhedra.polyhedron(v, Polyhedra.DefaultLibrary{Rational{Int64}}(GLPK.Optimizer))
            poly = ClutteredEnvPathOpt.gen_obstacle(points...)
            # display(poly)
            push!(obstacles, poly)

            obs_counter += 1
            points = Vector{Rational{Int}}[]
            if obs_counter == num_obs
                break
            end
        end
    end

    if !isempty(points)
        poly = ClutteredEnvPathOpt.gen_obstacle(points...)
        push!(obstacles, poly)
        obs_counter += 1
        points = Vector{Rational{Int}}[]
    end

    obstacles = ClutteredEnvPathOpt.gen_field(obstacles)   # maybe just apply remove_overlaps directly
    # Hardcoded
    if display_plot
        plot(title="Obstacles Seed $seed")
        ClutteredEnvPathOpt.plot_field(obstacles)   # requires an active plot (believe I have plot_field!, so need to edit plot_field)
        display(plot!())
        if save_plot
            # png("./test/obstacle files/Seed $seed")
            png("$(filename[1:end-3])")
        end
    end

    return obstacles
end

# Seeds by Sep 14
good_seeds = [100,99,98,97,96,95,94,93,92,89,86,83,80,78,77,75,73,70,69,68,67,66,64,62,58,57,56,55,54,53,51,50,48,47,46,45,43,42,39,37,36,35,34,33,32,30,26,25,24,23,22,20,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]
seeds = vcat(Vector(101:120), good_seeds, [201,202])  # all seeds

# Seeds by Sep 17, after choosing


seed = 107
num_obs = 1
file_name = "./test/obstacle files/Seed $seed.txt"

merge_faces = true;
partition = "CDT";
obstacles = obs_from_file(seed, num_obs, file_name)
obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(obstacles, "Obstacles Seed $seed Num Obs $num_obs", partition=partition, merge_faces=merge_faces)
skeleton = LabeledGraph(g)
all_faces = Set{Vector{Int64}}(union(obstacle_faces, free_faces))

# Quick analysis
plot(title="Obstacles Seed $seed");   # needs to be initiated
ClutteredEnvPathOpt.plot_field(obstacles)
# png("Obstacle Seed $seed Num Obs $num_obs")
display(plot!())
ClutteredEnvPathOpt.plot_points(points, vertices=skeleton.labels)
ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "Free Faces")
ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, plot_name = "Obstacle Faces", col = "dodgerblue")
ClutteredEnvPathOpt.plot_edges(skeleton, points, plot_name = "Skeleton", col="black", vertices=skeleton.labels)
feg_S = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
ClutteredEnvPathOpt.plot_edges(feg_S, points, plot_name = "Finite Element Graph of S", col="black", vertices=feg_S.labels)


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
savefig("Poster Obstacles Seed $seed Num Obs $num_obs.pdf")

# Plot faces
ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "Free Faces")
ClutteredEnvPathOpt.plot_faces(free_faces, points, plot_name = "Free Space Partition")
# plot!(title="",axis=([], false))
png("$partition Free Faces Seed $seed Num Obs $num_obs Merged Faces $merge_faces")
# png("$partition Free Faces Seed $seed Num Obs $num_obs Merged Faces $merge_faces No Axis")
savefig("Poster $partition Free Faces Seed $seed Num Obs $num_obs Merged Faces $merge_faces.pdf")

ClutteredEnvPathOpt.plot_faces(obstacle_faces, points, plot_name = "Obstacle Faces", col = "dodgerblue")
# plot!(title="",axis=([], false))
png("Obstacle Faces Seed $seed Num Obs $num_obs")

ClutteredEnvPathOpt.plot_faces(all_faces, points, plot_name = "All Faces", col = "red")
ClutteredEnvPathOpt.plot_faces(all_faces, points, plot_name = "All Faces", new_plot = false)
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

feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces))
ClutteredEnvPathOpt.plot_edges(feg, points, plot_name = "Finite Element Graph of S",col="black", vertices=feg_S.labels)
ClutteredEnvPathOpt.plot_edges(feg, points, plot_name = "G_S",col="black", vertices=feg_S.labels, with_labels=false)
# plot!(title="",axis=([], false))
png("$partition FEG of S Seed $seed Num Obs $num_obs Merged Faces $merge_faces")
# png("$partition FEG of S Seed $seed Num Obs $num_obs Merged Faces $merge_faces No Axis")
savefig("Poster $partition FEG of S Seed $seed Num Obs $num_obs Merged Faces $merge_faces.pdf")

# Plot conflict graph
conflict = LightGraphs.complement(feg.graph)
conflict_lg = LabeledGraph(conflict, feg.labels)
ClutteredEnvPathOpt.plot_edges(conflict_lg, points, plot_name = "Conflict Graph", col="black", with_labels=false)
png("$partition Conflict Graph Seed $seed Num Obs $num_obs Merged Faces $merge_faces")
savefig("Poster $partition Conflict Graph Seed $seed Num Obs $num_obs Merged Faces $merge_faces.pdf")

# feg_obs = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(obstacle_faces))
# ClutteredEnvPathOpt.plot_edges(feg_obs, points, plot_name = "Finite Element Graph Obstacle Faces")
# # plot!(title="",axis=([], false))
# png("$partition FEG Obstacles Seed $seed Num Obs $num_obs")

# feg_all = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces))
# ClutteredEnvPathOpt.plot_edges(feg_all, points, plot_name = "Finite Element Graph")
# # plot!(title="",axis=([], false))
# png("$partition FEG All Seed $seed Num Obs $num_obs")


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

# Plot biclique
ClutteredEnvPathOpt.plot_biclique_cover(feg, points, merged_cover; with_all=true, name="Poster Biclique")



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


#######################################################################################
########################### Biclique Cover Validity Tests #############################

# Biclique cover validity tester function
function biclique_cover_validity_tests(seed_range, num_obs_range; faces::String="free", with_plots::Bool=false, save_plots::Bool=false, partition="CDT", merge_faces::Bool=true, merge_BC::Bool=true)
    # fails = 0
    tuples = []
    for seed = seed_range
        for num_obs = num_obs_range
            println("On test Seed = $seed, Num_Obs = $num_obs")
            # Old with Random
            # obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(num_obs,"Obstacles Seed $seed Num Obs $num_obs",seed=seed, partition=partition, merge_faces=merge_faces)
            # skeleton = LabeledGraph(g)
            # all_faces = Set{Vector{Int64}}(union(obstacle_faces, free_faces))

            # New from file
            file_name = "./test/obstacle files/Seed $seed.txt"
            obstacles = obs_from_file(seed, num_obs, file_name, display_plot=false)
            obstacles, points, g, obstacle_faces, free_faces = ClutteredEnvPathOpt.plot_new(obstacles, "Obstacles Seed $seed Num Obs $num_obs", partition=partition, merge_faces=merge_faces)
            skeleton = LabeledGraph(g)
            all_faces = Set{Vector{Int64}}(union(obstacle_faces, free_faces))

            if faces == "all"
                cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, all_faces);
                feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(all_faces));
                merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
                if merge_BC == false
                    (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
                else
                    (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover);
                end
                # println("Difference: $(length(cover)-length(merged_cover)) | Merged: $(length(merged_cover)) | Original: $(length(cover))")
            else
                cover = ClutteredEnvPathOpt.find_biclique_cover(skeleton, free_faces);
                feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, ClutteredEnvPathOpt._find_face_pairs(free_faces));
                merged_cover = ClutteredEnvPathOpt.biclique_merger(cover, feg)
                if merge_BC == false
                    (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, cover);
                else
                    (valid_cov, size_diff, miss_edges, miss_edges_lg) = ClutteredEnvPathOpt._is_valid_biclique_cover_diff(feg, merged_cover);
                end
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
                        display(plot!())
                        if save_plots
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
    end
    return tuples
end

good_seeds = [100,99,98,97,96,95,94,93,92,89,86,83,80,78,77,75,73,70,69,68,67,66,64,62,58,57,56,55,54,53,51,50,48,47,46,45,43,42,39,37,36,35,34,33,32,30,26,25,24,23,22,20,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]
seeds = vcat(Vector(101:120), good_seeds, [201,202])  # all seeds

# TODO: Final seeds to test

seed_range = seeds   # 1:100
num_obs_range = 1:3
partition = "CDT"
merge_faces = false
merge_BC = false
# fail_tuples = biclique_cover_validity_tests(seed_range, num_obs_range, faces = "all", partition=partition)
fail_tuples = biclique_cover_validity_tests(seed_range, num_obs_range, faces = "free", with_plots=true, partition=partition, merge_faces=merge_faces, merge_BC=merge_BC)
# Note that whether the full cover or merged cover is checked is hardcoded in the function


######################################################################################
############################ Biclique Cover Merging Stats ############################

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
