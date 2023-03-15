using DataFrames
using LinearAlgebra
using Statistics
using Plots
using Pipe
using CSV
using StatsPlots
# using PyPlot
# using ClutteredEnvPathOpt
# using JuMP
# using Gurobi

####################################################################################
############################## Load individual files ###############################
# Load individual files
num_obs = 1
method = "merged"
partition = "CDT"
merge_face = false
file_name = "./Experiments/Solve Times/Solve Time Stats Seed Range All Num Obs $num_obs Method $method Partition $partition Merge Face $merge_face.txt"
# partition = "HP"
# file_name_hp = "Solve Time Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method $method Partition $partition.txt"

# # load via iostream
# raw_str = read(file_name, String)
# header = split(raw_str[3], '\t')
# header=[:seed, :num_obs, :solve_time, :rel_gap, :simplex_iter, :simplex_nodes, :num_vertices, :BC_merged, :BC_full, :num_free_faces, :num_free_face_ineq]
# io = IOBuffer(raw_str)
# df1 = CSV.File(io,
#                 delim='\t',
#                 ignorerepeated=true,
#                 header = header, # 3 for on line 3
#                 skipto = 4,
#                 ) |> 
#                 DataFrame

# Seed	Num_obs	Num_footsteps	Footsteps_used	Term_status	Obj_val	Solve_time	Rel_gap	Simplex_iterations	Barrier_iterations
# Nodes_explored	Num_vertices	BC_merged_size	BC_full_size	Num_free_faces	Num_free_face_inequalities	Last_f_cost	Between_f_cost	Trim_cost
header=[:seed, :num_obs, :num_footsteps, :footsteps_used, :term_status, :obj_val,
        :solve_time, :rel_gap, :simplex_iter, :barrier_iter, :simplex_nodes, :num_vertices,
        :BC_merged, :BC_full, :num_free_faces, :num_free_face_ineq, :last_f_cost, :between_f_cost, :trim_cost]
df = CSV.File(file_name,
                delim='\t',
                ignorerepeated=true,
                header = header, # 3 for on line 3
                skipto = 15,
                ) |> 
                DataFrame

##################################################################################
############################ Load files over loop ################################

ENV["COLUMNS"], ENV["LINES"] = 200, 15
header=[:seed, :num_obs, :num_footsteps, :footsteps_used, :term_status, :obj_val,
        :solve_time, :rel_gap, :simplex_iter, :barrier_iter, :simplex_nodes, :num_vertices,
        :BC_merged, :BC_full, :num_free_faces, :num_free_face_ineq, :last_f_cost, :between_f_cost, :trim_cost]

star_4 = union(Set([3,4,6,8,12,15,16,20,22,23,24,25,34,36,42,46,47,54,64,66,70,73,75,
        77,78,83,92,94,95,97,98,99,100]), Set(101:120))
star_3 = Set([1,5,9,13,14,17,33,37,39,45,51,53,55,62,80,89,96])
star_special = Set([201,202,203])
# seeds = union(star_4, star_3)
seeds = sort( collect( setdiff(union(star_4, star_3), Set(119) ) ) )

# seed_start = 1
# seed_end = 50
seed_range = seeds
num_obs_range = 1:3
partitions = ["CDT"]
# merge_faces = [true, false]
merge_faces = [false]
methods = ["merged", "full", "bigM"]
# methods = ["merged", "bigM"]
master_dict = Dict[]

for num_obs in num_obs_range
    dict = Dict()
    dictCDT = Dict()
    # dictHP = Dict()

    for partition in partitions
        # if partition == "CDT"
        # dictFM = Dict()
        dictFNM = Dict()

        for merge_face in merge_faces
            for method in methods
                file_name = "./Experiments/Solve Times/Solve Time Stats Seed Range All Num Obs $num_obs Method $method Partition $partition Merge Face $merge_face.txt"
                df = CSV.File(file_name,
                    delim='\t',
                    ignorerepeated=true,
                    header = header, # 3 for on line 3
                    skipto = 15,
                    ) |> 
                    DataFrame

                # CSV.write("DF $(file_name[1:end-4]).csv", df)
                
                dictMethod = Dict()
                dictMethod["filename"] = file_name
                dictMethod["df"] = df

                # if merge_face
                #     dictFM[method] = dictMethod
                # else
                dictFNM[method] = dictMethod
                # end
            end
        end

        # dictCDT["faces_merged"] = dictFM
        dictCDT["faces_not_merged"] = dictFNM
        dict["CDT"] = dictCDT
        # else
        #     if num_obs <= 2   # only have HP data in first two num_obs
        #         for method in methods
        #             file_name = "Solve Time Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method $method Partition $partition.txt"
        #             df = CSV.File(file_name,
        #                 delim='\t',
        #                 ignorerepeated=true,
        #                 header = header, # 3 for on line 3
        #                 skipto = 4,
        #                 ) |> 
        #                 DataFrame

        #             # CSV.write("DF $(file_name[1:end-4]).csv", df)

        #             dictMethod = Dict()
        #             dictMethod["filename"] = file_name
        #             dictMethod["df"] = df

        #             dictHP[method] = dictMethod
        #         end

        #         dict["HP"] = dictHP   
        #     end
        # end
    end

    # dict["CDT"] = dictCDT
    # if num_obs <= 2
    #     dict["HP"] = dictHP   # empty for num_obs 3 and 4
    # end
    push!(master_dict, dict)
end

num_obs = 1
master_dict[num_obs]["CDT"]["faces_not_merged"]["merged"]["df"]
master_dict[num_obs]["CDT"]["faces_not_merged"]["bigM"]["df"]
master_dict[num_obs]["CDT"]["faces_not_merged"]["full"]["df"]

############################################################################################
############################ Merge dictionaries one at a time ##############################
# # Adds "summary" DataFrame of the three methods

num_obs = 1
partition = "CDT"
faces = "faces_not_merged"
df = copy(master_dict[num_obs][partition][faces]["merged"]["df"])
rename!(df, :solve_time => :solve_time_merged, :term_status => :term_status_merged)
# Remove columns
# select!(df, Not([:rel_gap, :simplex_nodes, :simplex_iter, :num_obs]))
select!(df, :seed, :solve_time_merged, :term_status_merged, :num_vertices, :BC_merged, :BC_full, :num_free_faces, :num_free_face_ineq)
# Use : in rows to assure copy
df[:,:num_free_face_ineq] = master_dict[num_obs][partition][faces]["bigM"]["df"].num_free_face_ineq
# Add times from other methods
df[:, :solve_time_full] = master_dict[num_obs][partition][faces]["full"]["df"].solve_time
df[:, :term_status_full] = master_dict[num_obs][partition][faces]["full"]["df"].term_status
df[:, :solve_time_bigM] = master_dict[num_obs][partition][faces]["bigM"]["df"].solve_time
df[:, :term_status_bigM] = master_dict[num_obs][partition][faces]["bigM"]["df"].term_status
# Add biclique cover reduction stats
df.BC_reduction = map(i -> (df[i,:BC_full] - df[i,:BC_merged]) / df[i,:BC_full], 1:nrow(df))
# Reorder columns
# select!(df, :seed, :solve_time_merged, :solve_time_full, :solve_time_bigM,
#         Not([:num_free_faces, :num_free_face_ineq, :num_vertices, :BC_reduction]), :num_free_face_ineq,
#         :num_free_faces, :num_vertices, :BC_reduction)
select!(df, :seed, :solve_time_merged, :solve_time_bigM, :term_status_merged, :term_status_bigM,
        Not([:num_free_faces, :num_free_face_ineq, :num_vertices, :BC_reduction]), :num_free_face_ineq,
        :num_free_faces, :num_vertices, :BC_reduction)
# Add winner column
df.winner = map(i -> begin
            if df[i,:solve_time_merged] <= df[i,:solve_time_full] && df[i,:solve_time_merged] <= df[i,:solve_time_bigM] && (df[i, :term_status_merged] == "OPTIMAL")
                return "merged"
            elseif df[i,:solve_time_full] <= df[i,:solve_time_bigM] && (df[i, :term_status_full] == "OPTIMAL")
                return "full"
            # if (df[i,:solve_time_merged] < df[i,:solve_time_bigM]) && (df[i, :term_status_merged] == "OPTIMAL")
            #     return "merged"
            elseif (df[i, :term_status_bigM] == "OPTIMAL")
                return "bigM"
            else
                return "tie"
            end
        end,
        1:nrow(df))
# # Add winner column between merged and bigM
# df.winner_1v1 = map(i -> begin
#                 if df[i,:solve_time_merged] <= df[i,:solve_time_bigM]
#                     return "merged"
#                 else
#                     return "bigM"
#                 end
#             end,
#             1:nrow(df))
# Reorder columns
select!(df, :seed, :winner, :) # :winner_1v1, :)

# Store winner_count GroupedDataFrame
gdf = groupby(df, :winner)
cdf = combine(gdf, :winner => length => :winner_count)
master_dict[num_obs][partition][faces]["winner_count"] = cdf
# # Store winner_count_1v1 GroupedDataFrame
# gdf_1v1 = groupby(df, :winner_1v1)
# cdf_1v1 = combine(gdf_1v1, :winner_1v1 => length => :winner_count_1v1)
# master_dict[num_obs][partition][faces]["winner_count_1v1"] = cdf_1v1
# Store summary DataFrame
master_dict[num_obs][partition][faces]["summary"] = df

master_dict[num_obs][partition][faces]["summary"]
master_dict[num_obs][partition][faces]["winner_count"]
# master_dict[num_obs][partition][faces]["winner_count_1v1"]

############################################################################################
############################## Merge dictionaries over loop ################################
# Adds "summary" DataFrame of the three methods

for num_obs in 1:3
    for partition in ["CDT"]
        for faces in ["faces_not_merged"]
            df = copy(master_dict[num_obs][partition][faces]["merged"]["df"])
            rename!(df, :solve_time => :solve_time_merged, :term_status => :term_status_merged)
            # Remove columns
            # select!(df, Not([:rel_gap, :simplex_nodes, :simplex_iter, :num_obs]))
            select!(df, :seed, :solve_time_merged, :term_status_merged, :num_vertices, :BC_merged, :BC_full, :num_free_faces, :num_free_face_ineq)
            # Use : in rows to assure copy
            df[:,:num_free_face_ineq] = master_dict[num_obs][partition][faces]["bigM"]["df"].num_free_face_ineq
            # Add times from other methods
            df[:, :solve_time_full] = master_dict[num_obs][partition][faces]["full"]["df"].solve_time
            df[:, :term_status_full] = master_dict[num_obs][partition][faces]["full"]["df"].term_status
            df[:, :solve_time_bigM] = master_dict[num_obs][partition][faces]["bigM"]["df"].solve_time
            df[:, :term_status_bigM] = master_dict[num_obs][partition][faces]["bigM"]["df"].term_status
            # Add biclique cover reduction stats
            df.BC_reduction = map(i -> (df[i,:BC_full] - df[i,:BC_merged]) / df[i,:BC_full], 1:nrow(df))
            # Reorder columns
            # select!(df, :seed, :solve_time_merged, :solve_time_full, :solve_time_bigM,
            #         Not([:num_free_faces, :num_free_face_ineq, :num_vertices, :BC_reduction]), :num_free_face_ineq,
            #         :num_free_faces, :num_vertices, :BC_reduction)
            select!(df, :seed, :solve_time_merged, :solve_time_bigM, :term_status_merged, :term_status_bigM,
                    Not([:num_free_faces, :num_free_face_ineq, :num_vertices, :BC_reduction]), :num_free_face_ineq,
                    :num_free_faces, :num_vertices, :BC_reduction)
            # Add winner column
            df.winner = map(i -> begin
                        if df[i,:solve_time_merged] <= df[i,:solve_time_full] && df[i,:solve_time_merged] <= df[i,:solve_time_bigM] && (df[i, :term_status_merged] == "OPTIMAL")
                            return "merged"
                        elseif df[i,:solve_time_full] <= df[i,:solve_time_bigM] && (df[i, :term_status_full] == "OPTIMAL")
                            return "full"
                        # if (df[i,:solve_time_merged] < df[i,:solve_time_bigM]) && (df[i, :term_status_merged] == "OPTIMAL")
                        #     return "merged"
                        elseif (df[i, :term_status_bigM] == "OPTIMAL")
                            return "bigM"
                        else
                            return "tie"
                        end
                    end,
                    1:nrow(df))
            # # Add winner column between merged and bigM
            # df.winner_1v1 = map(i -> begin
            #                 if df[i,:solve_time_merged] <= df[i,:solve_time_bigM]
            #                     return "merged"
            #                 else
            #                     return "bigM"
            #                 end
            #             end,
            #             1:nrow(df))
            # Reorder columns
            select!(df, :seed, :winner, :) # :winner_1v1, :)

            # Store winner_count GroupedDataFrame
            gdf = groupby(df, :winner)
            cdf = combine(gdf, :winner => length => :winner_count)
            # @show cdf
            # println("\n")
            master_dict[num_obs][partition][faces]["winner_count"] = cdf
            # # Store winner_count_1v1 GroupedDataFrame
            # gdf_1v1 = groupby(df, :winner_1v1)
            # cdf_1v1 = combine(gdf_1v1, :winner_1v1 => length => :winner_count_1v1)
            # # @show cdf_1v1
            # # println("\n")
            # master_dict[num_obs][partition][faces]["winner_count_1v1"] = cdf_1v1
            # Store summary DataFrame
            master_dict[num_obs][partition][faces]["summary"] = df
        end
    end
end

num_obs = 3
master_dict[num_obs]["CDT"]["faces_not_merged"]["summary"]
master_dict[num_obs]["CDT"]["faces_not_merged"]["winner_count"]
# master_dict[num_obs]["CDT"]["faces_not_merged"]["winner_count_1v1"]


# Analysis on solve time gaps
for num_obs = 1:3
    df = master_dict[num_obs]["CDT"]["faces_not_merged"]["summary"]

    df_merged = df[df.winner .== "merged", :]
    # filter!(:term_status_bigM => ==("OPTIMAL"), df_merged)
    df_merged[:, :solve_time_bigM] - df_merged[:,:solve_time_merged]
    # ratios = (df_merged[:, :solve_time_bigM] - df_merged[:,:solve_time_merged]) ./ df_merged[:, :solve_time_bigM]
    println("\nNum Obs $num_obs Method merged")
    frac = df_merged[:,:solve_time_merged] ./ df_merged[:, :solve_time_bigM]
    display(frac)

    df_full = df[df.winner .== "full", :]
    # filter!(:term_status_bigM => ==("OPTIMAL"), df_full)
    df_full[:, :solve_time_bigM] - df_full[:,:solve_time_full]
    # ratios = (df_full[:, :solve_time_bigM] - df_full[:,:solve_time_full]) ./ df_full[:, :solve_time_bigM]
    println("\nNum Obs $num_obs Method full")
    frac = df_full[:,:solve_time_full] ./ df_full[:, :solve_time_bigM]
    display(frac)
end

df = master_dict[num_obs]["CDT"]["faces_not_merged"]["summary"]
df_merged = df[df.winner .== "merged", :]
# filter!(:term_status_bigM => ==("OPTIMAL"), df_merged)
df_merged[:, :solve_time_bigM] - df_merged[:,:solve_time_merged]

for num_obs = 1:3
    df_m = master_dict[num_obs]["CDT"]["faces_not_merged"]["merged"]["df"]
    df_f = master_dict[num_obs]["CDT"]["faces_not_merged"]["full"]["df"]
    df_b = master_dict[num_obs]["CDT"]["faces_not_merged"]["bigM"]["df"]

    temp = copy(df_m)
    select!(temp, :simplex_nodes)
    temp.simplex_nodes_full = df_f[:,:simplex_nodes]
    temp.simplex_nodes_bigM = df_b[:,:simplex_nodes]

    display(describe(temp))
    # combine(temp :simplex_nodes => mean, :simplex_nodes_full => mean, :simplex_nodes_bigM => mean)
end

############################################################################################
####################### Manually compute means and std for optimals ########################

master_dict[num_obs]["CDT"]["faces_not_merged"]
df = master_dict[num_obs]["CDT"]["faces_not_merged"]["merged"]["df"]

num_obs_range = 1:3
for num_obs in num_obs_range
    df = master_dict[num_obs]["CDT"]["faces_not_merged"]["merged"]["df"]
    # df_opt = df[df.term_status .== "OPTIMAL", :]
    gdf = groupby(df, :term_status)
    # opt_solve_mean = combine(gdf, :solve_time => mean)
    opt_solve_stats = combine(gdf, :solve_time => mean, :solve_time => x -> sqrt(var(x)), :solve_time => length)
    println("Num Obs $num_obs \n $opt_solve_stats")
end

for num_obs in num_obs_range
    df = master_dict[num_obs]["CDT"]["faces_not_merged"]["full"]["df"]
    # df_opt = df[df.term_status .== "OPTIMAL", :]
    gdf = groupby(df, :term_status)
    # opt_solve_mean = combine(gdf, :solve_time => mean)
    opt_solve_stats = combine(gdf, :solve_time => mean, :solve_time => x -> sqrt(var(x)), :solve_time => length)
    println("Num Obs $num_obs \n $opt_solve_stats")
end

for num_obs in num_obs_range
    df = master_dict[num_obs]["CDT"]["faces_not_merged"]["bigM"]["df"]
    # df_opt = df[df.term_status .== "OPTIMAL", :]
    gdf = groupby(df, :term_status)
    # opt_solve_mean = combine(gdf, :solve_time => mean)
    opt_solve_stats = combine(gdf, :solve_time => mean, :solve_time => x -> sqrt(var(x)), :solve_time => length)
    println("Num Obs $num_obs \n $opt_solve_stats")
end

##### BC cover sizes
num_obs = 3
df_sum = master_dict[num_obs]["CDT"]["faces_not_merged"]["summary"]
combine(df_sum, :BC_merged => mean, :BC_full => mean, :BC_reduction => mean, :num_free_faces => mean, :num_free_face_ineq => mean, :num_vertices => mean)
combine(df_sum, :BC_merged => x -> sqrt(var(x)), :BC_full => x -> sqrt(var(x)), :BC_reduction => x -> sqrt(var(x)),
        :num_free_faces => x -> sqrt(var(x)), :num_free_face_ineq => x -> sqrt(var(x)), :num_vertices => x -> sqrt(var(x)))

master_dict[num_obs]["CDT"]["faces_not_merged"]["winner_count"]

# file_name = "./Experiments/Problem Sizes/Problem Size Stats Seed Range All Seeds Num Obs $num_obs Method $method Partition $partition Merge Face $merge_face.txt"


##############################################################################################
###################### Load files with problem sizes into dictionaries #######################

ENV["COLUMNS"], ENV["LINES"] = 200, 15
header=[:seed, :num_obs, :disj_cont_var, :disj_bin_var, :disj_ineq_cons, :disj_eq_cons, :all_cont_var, :all_bin_var,
        :num_vertices, :BC_merged, :BC_full, :num_free_face_ineq, :num_free_faces]

star_4 = union(Set([3,4,6,8,12,15,16,20,22,23,24,25,34,36,42,46,47,54,64,66,70,73,75,
        77,78,83,92,94,95,97,98,99,100]), Set(101:120))
star_3 = Set([1,5,9,13,14,17,33,37,39,45,51,53,55,62,80,89,96])
star_special = Set([201,202,203])
# seeds = union(star_4, star_3)
seeds = sort( collect( setdiff(union(star_4, star_3), Set(119) ) ) )

seed_range = seeds
num_obs_range = 1:3
partitions = ["CDT"]
merge_faces = [false]
prob_methods = ["merged", "full", "bigM"]
master_sizes = Dict[]

for num_obs in num_obs_range
    dict = Dict()
    dictCDT = Dict()
    # dictHP = Dict()

    for partition in partitions
        # if partition == "CDT"
            # dictFM = Dict()
            dictFNM = Dict()

            for merge_face in merge_faces
                for method in prob_methods
                    file_name = "./Experiments/Problem Sizes/Problem Size Stats Seed Range All Seeds Num Obs $num_obs Method $method Partition $partition Merge Face $merge_face.txt"
                    df = CSV.File(file_name,
                        delim='\t',
                        ignorerepeated=true,
                        header = header, # 3 for on line 3
                        skipto = 3,
                        ) |> 
                        DataFrame

                    # CSV.write("DF $(file_name[1:end-4]).csv", df)
                    
                    dictMethod = Dict()
                    dictMethod["filename"] = file_name
                    dictMethod["df"] = df

                    # if merge_face
                    #     dictFM[method] = dictMethod
                    # else
                    dictFNM[method] = dictMethod
                    # end
                end
            end

            # dictCDT["faces_merged"] = dictFM
            dictCDT["faces_not_merged"] = dictFNM
            dict["CDT"] = dictCDT
        # else
        #     for method in prob_methods
        #         file_name = "Problem Size Stats Seed Range $seed_start to $seed_end Num Obs $num_obs Method $method Partition $partition.txt"
        #         df = CSV.File(file_name,
        #             delim='\t',
        #             ignorerepeated=true,
        #             header = header, # 3 for on line 3
        #             skipto = 4,
        #             ) |> 
        #             DataFrame

        #         # CSV.write("DF $(file_name[1:end-4]).csv", df)

        #         dictMethod = Dict()
        #         dictMethod["filename"] = file_name
        #         dictMethod["df"] = df

        #         dictHP[method] = dictMethod
        #     end

        #     dict["HP"] = dictHP   
        # end
    end

    push!(master_sizes, dict)
end

num_obs = 1
master_sizes[num_obs]["CDT"]["faces_not_merged"]["merged"]["df"]

############################################################################################
############################ Merge dictionaries one at a time ##############################
# # Adds "summary" DataFrame of the three methods

num_obs = 1
partition = "CDT"
faces = "faces_not_merged"
df = copy(master_sizes[num_obs][partition][faces]["merged"]["df"])
rename!(df, :disj_cont_var => :disj_cont_var_merged)
rename!(df, :disj_bin_var => :disj_bin_var_merged)
rename!(df, :disj_ineq_cons => :disj_ineq_cons_merged)
rename!(df, :disj_eq_cons => :disj_eq_cons_merged)
# Remove columns
select!(df, Not([:all_cont_var, :all_bin_var]))
# Use : in rows to assure copy
df[:,:num_free_face_ineq] = master_sizes[num_obs][partition][faces]["bigM"]["df"].num_free_face_ineq
# Add sizes from other methods
df[:, :disj_cont_var_full] = master_sizes[num_obs][partition][faces]["full"]["df"].disj_cont_var
df[:, :disj_bin_var_full] = master_sizes[num_obs][partition][faces]["full"]["df"].disj_bin_var
df[:, :disj_ineq_cons_full] = master_sizes[num_obs][partition][faces]["full"]["df"].disj_ineq_cons
df[:, :disj_eq_cons_full] = master_sizes[num_obs][partition][faces]["full"]["df"].disj_eq_cons

df[:, :disj_cont_var_bigM] = master_sizes[num_obs][partition][faces]["bigM"]["df"].disj_cont_var
df[:, :disj_bin_var_bigM] = master_sizes[num_obs][partition][faces]["bigM"]["df"].disj_bin_var
df[:, :disj_ineq_cons_bigM] = master_sizes[num_obs][partition][faces]["bigM"]["df"].disj_ineq_cons
df[:, :disj_eq_cons_bigM] = master_sizes[num_obs][partition][faces]["bigM"]["df"].disj_eq_cons
# Add biclique cover reduction stats
df.BC_reduction = map(i -> (df[i,:BC_full] - df[i,:BC_merged]) / df[i,:BC_full], 1:nrow(df))
# Reorder columns
select!(df, :seed, :disj_cont_var_merged, :disj_cont_var_full, :disj_cont_var_bigM,
        :disj_bin_var_merged, :disj_bin_var_full, :disj_bin_var_bigM,
        :disj_ineq_cons_merged, :disj_ineq_cons_full, :disj_ineq_cons_bigM,
        :disj_eq_cons_merged, :disj_eq_cons_full, :disj_eq_cons_bigM)
        # :num_vertices, :num_free_faces, :BC_merged, :BC_full,
        # :num_free_face_ineq, :BC_reduction)

############################################################################################
############################## Merge dictionaries over loop ################################
# Adds "summary" DataFrame of the three methods

for num_obs in 1:3
    for partition in ["CDT"]
        for faces in ["faces_not_merged"]
            df = copy(master_sizes[num_obs][partition][faces]["merged"]["df"])
            rename!(df, :disj_cont_var => :disj_cont_var_merged)
            rename!(df, :disj_bin_var => :disj_bin_var_merged)
            rename!(df, :disj_ineq_cons => :disj_ineq_cons_merged)
            rename!(df, :disj_eq_cons => :disj_eq_cons_merged)
            # Remove columns
            select!(df, Not([:all_cont_var, :all_bin_var]))
            # Use : in rows to assure copy
            df[:,:num_free_face_ineq] = master_sizes[num_obs][partition][faces]["bigM"]["df"].num_free_face_ineq
            # Add sizes from other methods
            df[:, :disj_cont_var_full] = master_sizes[num_obs][partition][faces]["full"]["df"].disj_cont_var
            df[:, :disj_bin_var_full] = master_sizes[num_obs][partition][faces]["full"]["df"].disj_bin_var
            df[:, :disj_ineq_cons_full] = master_sizes[num_obs][partition][faces]["full"]["df"].disj_ineq_cons
            df[:, :disj_eq_cons_full] = master_sizes[num_obs][partition][faces]["full"]["df"].disj_eq_cons

            df[:, :disj_cont_var_bigM] = master_sizes[num_obs][partition][faces]["bigM"]["df"].disj_cont_var
            df[:, :disj_bin_var_bigM] = master_sizes[num_obs][partition][faces]["bigM"]["df"].disj_bin_var
            df[:, :disj_ineq_cons_bigM] = master_sizes[num_obs][partition][faces]["bigM"]["df"].disj_ineq_cons
            df[:, :disj_eq_cons_bigM] = master_sizes[num_obs][partition][faces]["bigM"]["df"].disj_eq_cons
            # Add biclique cover reduction stats
            df.BC_reduction = map(i -> (df[i,:BC_full] - df[i,:BC_merged]) / df[i,:BC_full], 1:nrow(df))
            # Reorder columns
            select!(df, :seed, :disj_cont_var_merged, :disj_cont_var_full, :disj_cont_var_bigM,
                    :disj_bin_var_merged, :disj_bin_var_full, :disj_bin_var_bigM,
                    :disj_ineq_cons_merged, :disj_ineq_cons_full, :disj_ineq_cons_bigM,
                    :disj_eq_cons_merged, :disj_eq_cons_full, :disj_eq_cons_bigM)
                    # :num_vertices, :num_free_faces, :BC_merged, :BC_full,
                    # :num_free_face_ineq, :BC_reduction)

            # Store summary DataFrame
            master_sizes[num_obs][partition][faces]["summary"] = df
        end
    end
end

############################################################################################
########################### Manually compute problem size means ############################

num_obs = 3
df = master_sizes[num_obs]["CDT"]["faces_not_merged"]["summary"]
temp = combine(df, :disj_cont_var_merged => mean, :disj_cont_var_full => mean, :disj_cont_var_bigM => mean,
        :disj_bin_var_merged => mean, :disj_bin_var_full => mean, :disj_bin_var_bigM => mean,
        :disj_ineq_cons_merged => mean, :disj_ineq_cons_full => mean, :disj_ineq_cons_bigM => mean,
        :disj_eq_cons_merged => mean, :disj_eq_cons_full => mean, :disj_eq_cons_bigM => mean)
temp[:,1:3]
temp[:,4:6]
temp[:,7:9]
temp[:,10:end]