
# Obstacle Selection
# All generated with num_obs = 3

using ClutteredEnvPathOpt
using Plots

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

# Original list from Random
# one_obs_seeds = Set([24,30,33,42])
# two_obs_seeds = Set([4,6,8,12,14,15,19,20,21,22,25,28,34,36,37,39,41,46,47,49,51,53,63,64,65,72,73,75,76,77,79,81,82,83,84,85,87,88,90,91,94,98,99])
# three_obs_seeds = setdiff( setdiff( Set(1:100), two_obs_seeds), one_obs_seeds)

# Sep 14
good_seeds = Set([100,99,98,97,96,95,94,93,92,89,86,83,80,78,77,75,73,70,69,68,67,66,64,62,58,57,56,55,54,53,51,50,
                    48,47,46,45,43,42,39,37,36,35,34,33,32,30,26,25,24,23,22,20,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1])
good_seeds_vec = [100,99,98,97,96,95,94,93,92,89,86,83,80,78,77,75,73,70,69,68,67,66,64,62,58,57,56,55,54,53,51,50,
                    48,47,46,45,43,42,39,37,36,35,34,33,32,30,26,25,24,23,22,20,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]
removed_seeds = Set([19,21,63,65,90,91])   # were originally in good_seeds
# Up-to-date notes
# 3 and 4 can see effect of placement of an obstacle (that's not in the shortest path)
# 3 and 8 can see effect of size of an obstacle (that's not in the shortest path)
# 25 and 99 see effect of changing a pentagon to a square
# 22 and 97 can see effect of the same obstacle being on boundary vs not on boundary
# 22 abd 95 can see effect of removing a wide angle vertex with wide angle
# 9 and 13 complement each other
# good_seeds that became another seed: 
    # 19->104,21->109 and 21->110,63->108,65->106; 90->201,91->202

# Original list from Random
meh_seeds = Set([88,87,84,82,76,72,71,60,59,52,49,44,41,38,28,27])
bad_seeds = Set([85,81,79,74,61,40,31,29])
# Files for these seeds have been deleted

Set(1:100) == union(good_seeds, meh_seeds, bad_seeds, removed_seeds)

#############################################################################################
#################################### Some notes on seeds ####################################

# Seeds with boundary obstacles
union(Set([16,77,97,115,116,118,119,120]), Set(101:111))

# Seeds with wide angles
Set([15,22,26,30,39,51,53,56,64,95,108])  # Remember that resulting polytopes must still have 4+ vertices

# Seeds good for seeing how small close obtacles compare vs their convex hull
# Maybe do 5 to 10 cases of the above to compare (name them in the 200 range)
Set([12,13,15,20,22,35,66,73,83,93,94,96])  # 22,35,66,93,96 combine two obs;

# Seeds with lots of vertices (or one obstacle with 6+)
Set([6,16,24,30,33,36,46,51,64,73,75])   #  73, 75

# Seeds with small obstacles
Set([22,95,97,112,113,114,115,117])

# Huge obstacles or convex hull cases
Set([201,202])

# Pure square obstacles
Set([107,111,115,116,117,119])

# From top
all_seeds = vcat(Vector(101:120), good_seeds_vec, [201,202])

# This is the final set of seeds. We will use star_4 and star_3. 
star_4 = union(Set([3,4,6,8,12,15,16,20,22,23,24,25,34,36,42,46,47,54,64,66,70,73,75,
                77,78,83,92,94,95,97,98,99,100]), Set(101:120))
star_3 = Set([1,5,9,13,14,17,33,37,39,45,51,53,55,62,80,89,96])
star_2 = Set([2,7,10,11,18,30,35,48,56,57,58,67,68,69,86,93])
star_1 = Set([26,32,43,50])
star_special = Set([201,202])

# Check if seeds match
all_seeds = union(star_4, star_3, star_2, star_1, star_special)
all_seeds_2 = vcat(Vector(101:120), good_seeds_vec, [201,202])   # from top
Set(all_seeds_2) == all_seeds


#############################################################################################
################################ Display obstacles from files ###############################


num_obs = 3
seeds = sort!( collect( union(star_4, star_3)))
seeds = [118]

for seed in seeds
    file_name = "./test/obstacle files/Seed $seed.txt"
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

    # Remember to run gen_field or remove_overlaps on the obstacle vector to remove overlaps
    obstacles = ClutteredEnvPathOpt.gen_field(obstacles)
    plot(title="Obstacles Seed $seed")
    ClutteredEnvPathOpt.plot_field(obstacles)   # requires an active plot (believe I have plot_field!, so need to edit plot_field)
    display(plot!())
    # png("./test/obstacle files/Seed $seed")
end


# Plot grid points
plot(legend=false)
x = [i//21 for i = 0:21]
for j = 0:21
    y = [j//21 for i = 0:21]
    scatter!(x,y)
end
display(plot!())


# TODO: Write code for reading Floats or Ints


#############################################################################################
################################# Some previous seed notes ##################################

# Obstacles close to the starting point
# Seed    Description
# 1   Two close together, one behind them a bit back. Should walk above both
# 16

# Obstacles "on" shortest path
# Seed    Description
# 2   Path should make an arc (concave) in between the Obstacles
# 3   Similar to seed 2 but it bit rotated clockwise
# 10
# 11
# 15  Two obstacles. One is large with 5 vertices
# 15  Two obstacles. One is large with 6 vertices
# 16
# 18

# Obstacles near the goal position
# Seed    Description
# 4   Two obstacle near the goal. They are really close together. One has 5 vertices
# 11

# Obstacles are spaced out
# Seed    Description
# 5   Obstacles make a capital gamma shape, with one being close to goal position
# 6   Two obstacles on the lower portion of the cell. One has 6 vertices
# 18

# Robot should walk under
# Seed    Description
# 7   Can also walk through
# 12  Two ostacles. One has 5 vertices
# 13

# Robot should walk below
# Seed    Description
# 9   

# Obstacles not really in the way
# 9   
# 13  
# 14  Should walk between
# 19  Two obstacles. One has 5 vertices

# Obstacles are at top and at bottom, or left and right, or robot should walk between
# Seed    Description
# 6
# 8   Two obstacles
# 18
# 19
# 20
