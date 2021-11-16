"""
    find_biclique_cover(skeleton, faces)

Given a skeleton of a planar graph and its faces, finds a biclique cover for the complement
of the respective finite element graph with a divide-and-conquer approach: construct a biclique 
from a separator for the current finite element graph and call the function recursively on the 
induced subskeletons and subfaces of A ∪ C and B ∪ C.
"""
function find_biclique_cover(skeleton::LabeledGraph{T}, faces::Set{Vector{T}})::Set{Pair{Set{T}, Set{T}}} where {T}
    face_pairs = _find_face_pairs(faces)
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)

    if ClutteredEnvPathOpt.lg_is_complete(feg)
        return Set{Pair{Set{T}, Set{T}}}()
    end

    (C, A, B) = _find_feg_separator_lt_no_empty(skeleton, face_pairs)
    # println("Set A: $A\nSet B: $B\nSet C: $C")

    if !isempty(A) && !isempty(B)
        node = Set([A => B])
    else
        node = Set{Pair{Set{T}, Set{T}}}()
        println("A or B are empty.\nSkeleton vertices are $(skeleton.labels)")
    end

    if !isempty(union(A,C))
        skeleton_ac, faces_ac = _find_skeleton_faces(union(A, C), skeleton, faces)
        if !isempty(faces_ac)
            left = find_biclique_cover(skeleton_ac, faces_ac)
        else
            left = Set{Pair{Set{T}, Set{T}}}()
        end
    else
        left = Set{Pair{Set{T}, Set{T}}}()
    end

    if !isempty(union(B,C))
        skeleton_bc, faces_bc = _find_skeleton_faces(union(B, C), skeleton, faces)
        if !isempty(faces_bc)
            right = find_biclique_cover(skeleton_bc, faces_bc)
        else
            right = Set{Pair{Set{T}, Set{T}}}()
        end
    else
        right = Set{Pair{Set{T}, Set{T}}}()
    end

    return union(node, left, right)
end

"""
    find_biclique_cover_debug(skeleton, faces, points, obstacles, file_name)

Same as find_biclique_cover, with the added features of plotting the new skeletons and
recording separator information at each call.
"""
function find_biclique_cover_debug(skeleton::LabeledGraph{T}, faces::Set{Vector{T}}, points::Vector{Any}, obstacles, file_name)::Set{Pair{Set{T}, Set{T}}} where {T}
    ClutteredEnvPathOpt.plot_edges(skeleton, points, obstacles)
    f = open(file_name, "a")
    for (i,face) in enumerate(faces)
        write(f, "Face $i: $face\n")
        # println("Face $i: $face")
    end
    flush(f)
    close(f)
    
    face_pairs = _find_face_pairs(faces)
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)

    if ClutteredEnvPathOpt.lg_is_complete(feg)
        f = open(file_name, "a")
        write(f, "FEG is complete.\n")
        flush(f)
        close(f)

        # println("FEG is complete.")
        return Set{Pair{Set{T}, Set{T}}}()
    end

    (C, A, B) = _find_feg_separator_lt_no_empty(skeleton, face_pairs)

    if !isempty(A) && !isempty(B)
        node = Set([A => B])
    else
        node = Set{Pair{Set{T}, Set{T}}}()
        f = open(file_name, "a")
        write(f, "A or B is empty.\n")
        flush(f)
        close(f)

        # println("A or B is empty.")
    end

    if !isempty(union(A,C))
        skeleton_ac, faces_ac = _find_skeleton_faces(union(A, C), skeleton, faces)
        if !isempty(faces_ac)  # I think this will always be nonempty of A ∪ C is nonempty
            f = open(file_name, "a")
            write(f, "\nPlot \nSet A: $A\nSet B: $B\nSet C: $C\nA ∪ C\n")
            flush(f)
            close(f)

            # println("\nSet A: $A\nSet B: $B\nSet C: $C\nA ∪ C")
            left = find_biclique_cover_debug(skeleton_ac, faces_ac, points, obstacles, file_name)
        else
            left = Set{Pair{Set{T}, Set{T}}}()
            f = open(file_name, "a")
            write(f, "faces_ac is empty.\n")
            flush(f)
            close(f)

            # println("faces_ac is empty.")
        end
    else
        left = Set{Pair{Set{T}, Set{T}}}()
        f = open(file_name, "a")
        write(f, "A ∪ C is empty.\n")
        flush(f)
        close(f)

        # println("A ∪ C is empty.")
    end

    if !isempty(union(B,C))
        skeleton_bc, faces_bc = _find_skeleton_faces(union(B, C), skeleton, faces)
        if !isempty(faces_bc)
            f = open(file_name, "a")
            write(f, "\nPlot \nSet A: $A\nSet B: $B\nSet C: $C\nB ∪ C\n")
            flush(f)
            close(f)

            # println("\nSet A: $A\nSet B: $B\nSet C: $C\nB ∪ C")
            right = find_biclique_cover_debug(skeleton_bc, faces_bc, points, obstacles, file_name)
        else
            right = Set{Pair{Set{T}, Set{T}}}()
            f = open(file_name, "a")
            write(f, "faces_bc is empty.\n")
            flush(f)
            close(f)

            # println("faces_bc is empty.")
        end
    else
        right = Set{Pair{Set{T}, Set{T}}}()
        f = open(file_name, "a")
        write(f, "B ∪ C is empty.\n")
        flush(f)
        close(f)

        # println("B ∪ C is empty.")
    end

    return union(node, left, right)
end

function find_biclique_cover_one_iter(skeleton::LabeledGraph{T}, faces::Set{Vector{T}}) where {T}
    face_pairs = _find_face_pairs(faces)
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)

    if ClutteredEnvPathOpt.lg_is_complete(feg)
        return Set{Pair{Set{T}, Set{T}}}()
    end

    (C, A, B) = _find_feg_separator_lt_no_empty(skeleton, face_pairs)
    # println("Set A: $A\nSet B: $B\nSet C: $C")

    # if !isempty(A) && !isempty(B)
    #     node = Set([A => B])
    # else
    #     node = Set{Pair{Set{T}, Set{T}}}()
    #     println("A or B are empty.\nSkeleton vertices are $(skeleton.labels)")
    # end

    if !isempty(union(A,C))
        skeleton_ac, faces_ac = _find_skeleton_faces(union(A, C), skeleton, faces)
    #     if !isempty(faces_ac)
    #         left = find_biclique_cover(skeleton_ac, faces_ac)
    #     else
    #         left = Set{Pair{Set{T}, Set{T}}}()
    #     end
    # else
    #     left = Set{Pair{Set{T}, Set{T}}}()
    end

    if !isempty(union(B,C))
        skeleton_bc, faces_bc = _find_skeleton_faces(union(B, C), skeleton, faces)
    #     if !isempty(faces_bc)
    #         right = find_biclique_cover(skeleton_bc, faces_bc)
    #     else
    #         right = Set{Pair{Set{T}, Set{T}}}()
    #     end
    # else
    #     right = Set{Pair{Set{T}, Set{T}}}()
    end

    # return union(node, left, right)
    return (C,A,B), skeleton_ac, faces_ac, skeleton_bc, faces_bc
end

"""
    biclique_merger(cover, feg)

Given a biclique cover, merges bicliques together when possible for a smaller cover.
"""
function biclique_merger(cover::Set{Pair{Set{T}, Set{T}}}, feg::LabeledGraph{T})::Set{Pair{Set{T},Set{T}}} where {T}
    # Sort bicliques by size so that we attempt merging smaller bicliques first
    cover = sort(collect(cover), by=(s->length(s.first)+length(s.second)))
    
    # Loop until we are unable to find any merges
    merge_found = true
    while merge_found
        for i = 1:length(cover)
            for j = i+1:length(cover)
                # TODO: Check if this is right
                # if (!isempty(intersect(cover[i].first, cover[j].first)) && !isempty(intersect(cover[i].first, cover[j].second))) ||
                #     (!isempty(intersect(cover[i].second, cover[j].first)) && !isempty(intersect(cover[i].second, cover[j].second)))
                #     continue
                # end


                potential_merge_found = false
                merge_found = false
                A1_adj_A2 = false
                A1_adj_B2 = false
                B1_adj_A2 = false
                B1_adj_B2 = false
                for u in cover[i].first
                    for v in cover[j].first
                        if LightGraphs.has_edge(feg.graph, feg.labels[u], feg.labels[v])
                            A1_adj_A2 = true
                            break
                        end
                    end
                    if A1_adj_A2
                        break
                    end
                end

                for u in cover[i].first
                    for v in cover[j].second
                        if LightGraphs.has_edge(feg.graph, feg.labels[u], feg.labels[v])
                            A1_adj_B2 = true
                            break
                        end
                    end
                    if A1_adj_B2
                        break
                    end
                end

                if A1_adj_A2 && A1_adj_B2
                    continue  # cannot merge the bicliques in this case
                end

                if !A1_adj_A2
                    for u in cover[i].second 
                        for v in cover[j].second
                            if LightGraphs.has_edge(feg.graph, feg.labels[u], feg.labels[v])
                                B1_adj_B2 = true
                                break
                            end
                        end
                        if B1_adj_B2
                            break
                        end
                    end
                end
                
                if !A1_adj_B2
                    for u in cover[i].second 
                        for v in cover[j].first
                            if LightGraphs.has_edge(feg.graph, feg.labels[u], feg.labels[v])
                                B1_adj_A2 = true
                                break
                            end
                        end
                        if B1_adj_A2
                            break
                        end
                    end
                end

                # Two possible merges. If both are possible, will need to see if it matters which merge we choose
                A = Set([])
                B = Set([])
                if !A1_adj_A2 && !B1_adj_B2
                    # A = A1 ∪ B2
                    A = union(cover[i].first, cover[j].second)
                    # B = B1 ∪ A2
                    B = union(cover[i].second, cover[j].first)

                    potential_merge_found = true
                elseif !A1_adj_B2 && !B1_adj_A2
                    # A = A1 ∪ A2
                    A = union(cover[i].first, cover[j].first)
                    # B = B1 ∪ B2
                    B = union(cover[i].second, cover[j].second)
                    
                    potential_merge_found = true
                end

                # Second condition is to prevent vertices from being in both A and B
                if potential_merge_found && isempty(intersect(A,B)) # don't need second condition if above is correct
                    merge_found = true
                    deleteat!(cover, j)
                    deleteat!(cover, i)
                    push!(cover, A => B)
                    sort!(cover, by=(s->length(s.first)+length(s.second)))
                    break
                end
            
            end
            if merge_found
                break
            end
        end

    end

    return Set(cover)
end

struct BicliqueCoverTree
    biclique
    left
    right
end

function find_biclique_cover_as_tree(skeleton::LabeledGraph{T}, faces::Set{Vector{T}}) where {T}
    face_pairs = ClutteredEnvPathOpt._find_face_pairs(faces)
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)

    if ClutteredEnvPathOpt.lg_is_complete(feg)
        return nothing
    end

    (C, A, B) = ClutteredEnvPathOpt._find_feg_separator_lt_no_empty(skeleton, face_pairs)

    if !isempty(A) && !isempty(B)
        node = A => B
    else
        node = nothing
    end
    if !isempty(union(A,C))
        skeleton_ac, faces_ac = ClutteredEnvPathOpt._find_skeleton_faces(union(A, C), skeleton, faces)
        if !isempty(faces_ac)
            left = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton_ac, faces_ac)
        else
            left = nothing
        end
    else
        left = nothing
    end
    if !isempty(union(B,C))
        skeleton_bc, faces_bc = ClutteredEnvPathOpt._find_skeleton_faces(union(B, C), skeleton, faces)
        if !isempty(faces_bc)
            right = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton_bc, faces_bc)
        else
            right = nothing
        end
    else
        right = nothing
    end

    return BicliqueCoverTree(node, left, right)
end

function flatten(tree)
    if (tree === nothing)
        return []
    end

    return union([tree.biclique], flatten(tree.left), flatten(tree.right))
end

function getEdges(tree, cover)
    res = ""

    if (tree === nothing)
        return res
    end

    node_num = findfirst(x -> x === tree.biclique, cover)

    if (tree.left !== nothing) 
        left_num = findfirst(x -> x === tree.left.biclique, cover)
        res = string(res, node_num, " -> ", left_num, " ", getEdges(tree.left, cover))
    end

    if (tree.right !== nothing) 
        right_num = findfirst(x -> x === tree.right.biclique, cover)
        res = string(res, node_num, " -> ", right_num, " ", getEdges(tree.right, cover))
    end

    return res
end

function tree2digraph(tree)
    cover = collect(flatten(tree))

    edges = getEdges(tree, cover)

    return (cover, string("Digraph G { ", edges, " }"))
end

"""
    _find_skeleton_faces(vertices, old_skeleton, old_faces)

Find the skeleton and faces of a finite element subgraphgraph given a subset of
vertices, a finite element graph's skeleton and the faces of its planar
embedding. Returns a (skeleton, face vector) tuple.

Note: ensure the faces are passed in as a set of clockwise-ordered vertices
"""
function _find_skeleton_faces(vertices::Set{T}, old_skeleton::LabeledGraph{T}, old_faces::Set{Vector{T}})::Tuple{LabeledGraph{T}, Set{Vector{T}}} where {T}
    skeleton = copy(old_skeleton)
    for vertex in keys(old_skeleton.labels)
        if !in(vertex, vertices)
            ClutteredEnvPathOpt.lg_rem_vertex!(skeleton, vertex)
        end
    end

    # Remove obsolete vertices from each face and return faces with at least 1 vertex
    faces = @pipe map(face -> filter(vertex -> in(vertex, vertices), face), collect(old_faces)) |>
                filter(face -> length(face) > 0, _)
    
    # If a new face is a subset of another new face, don't include it
    new_faces = []
    removed_indices = []
    for i = 1:length(faces)
        indicator = 0
        for j = 1:length(faces)
            if j == i || j in removed_indices
                continue
            end
            if intersect(Set(faces[i]), Set(faces[j])) == Set(faces[i]) && length(faces[j]) >= length(faces[i])
                indicator += 1
                break
            end
        end
        if indicator != 0
            push!(removed_indices, i)
        else
            push!(new_faces, faces[i])
        end
    end

    # Close the new subface
    for face in new_faces
        if length(face) > 1
            for i = 2:length(face)
                ClutteredEnvPathOpt.lg_add_edge!(skeleton, face[i-1], face[i])
            end
            ClutteredEnvPathOpt.lg_add_edge!(skeleton, face[end], face[1])
        end
    end

    return skeleton, Set(new_faces)
end

"""
    _find_face_pairs(faces)

Convert a set of list of vertices included in a face to a set of sets of edges
comprising a face.
"""
function _find_face_pairs(faces::Set{Vector{T}})::Set{Set{Pair{T, T}}} where {T}
    face_pairs = Set{Set{Pair{T, T}}}()
    for face in faces
        pairs = Set{Pair{T, T}}()

        for i in 2:length(face)
            push!(pairs, face[i - 1] => face[i])
        end

        push!(pairs, face[end] => face[1])

        push!(face_pairs, pairs)
    end

    return face_pairs
end

"""
    _find_feg_separator_lt_no_empty(skeleton, face_pairs)

Find the separator of a finite element graph, repeating the process if for some
given root either A or B is empty. Returns the separator for which |A|*|B| is largest.
"""
function _find_feg_separator_lt_no_empty(skeleton::LabeledGraph{T}, face_pairs::Set{Set{Pair{T, T}}})::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)
    # for root in keys(skeleton.labels)
    #     (C, A, B) = find_feg_separator_lt(skeleton, face_pairs, root)
    #     (C, A, B) = pp_expell(feg, C, A, B)

    #     if !isempty(A) && !isempty(B)
    #         return (C, A, B)
    #     end
    # end

    # return (Set([]), Set([]), Set([]))
    #####################################################################################################

    separators = @pipe map(root -> find_feg_separator_lt(skeleton, face_pairs, root), collect(keys(skeleton.labels))) |>
        map(sep -> pp_expell(feg, sep[1], sep[2], sep[3]), _)

    balances = map(separator -> begin
        if isempty(separator[2]) || isempty(separator[3]) # A or B empty is no bueno
            return -1.0
        elseif isempty(separator[1]) # C empty, A and B are not, prioritize this case
            return 1.0 + min(length(separator[2]) / length(separator[3]), length(separator[3]) / length(separator[2]))
        else # None are empty
            return 1.0  # the separator that finds the most edges will be chosen
            # return min(length(separator[2]) / length(separator[3]), length(separator[3]) / length(separator[2]))
        end
    end, separators)

    max_balance = max(balances...)
    # println("Max Balance: $max_balance")

    if max_balance < -0.5
        # If at least one of A or B is empty, the subskeleton would be the same and we'll have infinite recursion
        # In such a case, our SEPARATOR ALGORITHM and PP_EXPELL cannot produce two nonempty A and B
        # Return empty sets since there'll be no biclique
        return (Set([]), Set([]), Set([]))
    else
        # Return the biclique that makes the most edges (e.g. |A|*|B| is largest)
        select = filter(b -> balances[b] == max_balance, 1:length(balances))
        max_edges = 0
        best_index = 0
        for i in select
            if length(separators[i][2])*length(separators[i][3]) > max_edges
                max_edges = length(separators[i][2])*length(separators[i][3])
                best_index = i
            end
        end
        return separators[best_index]
    end
    
end

"""
    _is_valid_biclique_cover(lg, cover)

Tests whether or not a cover is a valid biclique cover of a given graph.
"""
function _is_valid_biclique_cover(lg::LabeledGraph{T}, cover::Set{Pair{Set{T}, Set{T}}})::Bool where {T}
    e_bar = LightGraphs.edges(LightGraphs.complement(lg.graph))

    edges = @pipe map(pair -> _cartesian_product(pair.first, pair.second), collect(cover)) |> reduce(union!, _, init=Set{Pair{T, T}}())
    test = LightGraphs.SimpleGraph(LightGraphs.nv(lg.graph))

    for edge in edges
        # edge.first/second will be our nodes, so we map them to the nodes they are in the graph
        LightGraphs.add_edge!(test, lg.labels[edge.first], lg.labels[edge.second])
    end

    return LightGraphs.edges(test) == e_bar
end

"""
    _is_valid_biclique_cover_diff(lg, cover)

Tests whether or not a cover is a valid biclique cover of a given graph. Returns 
their difference in cardinality (positive => missing edges, negative => surplus of edges)
and a LabeledGraph containing any missing edges.
"""
function _is_valid_biclique_cover_diff(lg::LabeledGraph{T}, cover::Set{Pair{Set{T}, Set{T}}})::Tuple{Bool,Number,Set{LightGraphs.SimpleGraphs.SimpleEdge{T}},LabeledGraph{T}} where {T}
    e_bar = LightGraphs.edges(LightGraphs.complement(lg.graph))

    edges = @pipe map(pair -> _cartesian_product(pair.first, pair.second), collect(cover)) |> reduce(union!, _, init=Set{Pair{T, T}}())
    test = LightGraphs.SimpleGraph(LightGraphs.nv(lg.graph))
    rev = _reverse_labels(lg.labels)

    for edge in edges
        # edge.first/second will be our nodes, so we map them to the nodes they are in the graph
        LightGraphs.add_edge!(test, lg.labels[edge.first], lg.labels[edge.second])
    end

    # Create a LabeledGraph containing only the missing edges
    set_diff = setdiff(Set(e_bar), Set(LightGraphs.edges(test)))
    missing_e_graph = LightGraphs.SimpleGraph(LightGraphs.nv(lg.graph))
    missing_e_lg = LabeledGraph{Int}(missing_e_graph, lg.labels)
    for edge in set_diff
        lg_add_edge!(missing_e_lg, rev[edge.src], rev[edge.dst])
    end

    return LightGraphs.edges(test) == e_bar, length(e_bar) - LightGraphs.ne(test), set_diff, missing_e_lg
end

"""
    _cartesian_product(A, B)

Returns the cartesian product of two sets
"""
function _cartesian_product(A::Set{V}, B::Set{W})::Set{Pair{V, W}} where {V, W}
    res = Set{Pair{V, W}}()

    for a in A
        for b in B
            push!(res, Pair(a, b))
        end
    end

    return res
end

"""
    _find_feg_separator_lt_no_empty_visual(skeleton, face_pairs, locations)
Finds (C, A, B) and plots the edges found between A and B.
"""
function _find_feg_separator_lt_no_empty_visual(skeleton::LabeledGraph{T}, face_pairs::Set{Set{Pair{T, T}}}, points::Vector{Any})::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    (C,A,B) = _find_feg_separator_lt_no_empty(skeleton, face_pairs)

    if isempty(A) || isempty(B)
        return (C,A,B)
    end

    # Create edges between A and B
    new_edges = ClutteredEnvPathOpt._cartesian_product(A, B)

    Plots.plot()
    Plots.scatter!(map(point -> point.first, points), map(point -> point.second, points))

    for a in A
        Plots.scatter!([points[a].first], [points[a].second], color=:red)
    end
    for b in B
        Plots.scatter!([points[b].first], [points[b].second], color=:yellow)
    end

    for edge in new_edges
        (x1,y1) = points[edge.first]
        (x2,y2) = points[edge.second]
        Plots.plot!([x1, x2], [y1, y2], color=:green, linestyle=:dash)
    end

    display(Plots.plot!(legend=:false,xlims=(-.05,1.05),ylims=(-0.05,1.05), title="Biclique"))

    return (C,A,B)
end

"""
    _find_biclique_cover_visual(skeleton, faces, locations)
Finds a biclique cover and plots the edges of each biclique.
"""
function _find_biclique_cover_visual(skeleton::LabeledGraph{T}, faces::Set{Vector{T}}, points::Vector{Any})::Set{Pair{Set{T}, Set{T}}} where {T}
    face_pairs = _find_face_pairs(faces)
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)

    if ClutteredEnvPathOpt.lg_is_complete(feg)
        return Set{Pair{Set{T}, Set{T}}}()
    end

    (C, A, B) = _find_feg_separator_lt_no_empty_visual(skeleton, face_pairs, points)
    # println("Set A: $A\nSet B: $B\nSet C: $C")

    if !isempty(A) && !isempty(B)
        node = Set([A => B])
    else
        node = Set{Pair{Set{T}, Set{T}}}()
    end
    if !isempty(union(A,C))
        skeleton_ac, faces_ac = _find_skeleton_faces(union(A, C), skeleton, faces)
        if !isempty(faces_ac)
            # println("\nSet A: $A\nSet B: $B\nSet C: $C\nA ∪ C")
            left = _find_biclique_cover_visual(skeleton_ac, faces_ac, points)
        else
            left = Set{Pair{Set{T}, Set{T}}}()
        end
    else
        left = Set{Pair{Set{T}, Set{T}}}()
    end
    if !isempty(union(B,C))
        skeleton_bc, faces_bc = _find_skeleton_faces(union(B, C), skeleton, faces)
        if !isempty(faces_bc)
            # println("\nSet A: $A\nSet B: $B\nSet C: $C\nB ∪ C")
            right = _find_biclique_cover_visual(skeleton_bc, faces_bc, points)
        else
            right = Set{Pair{Set{T}, Set{T}}}()
        end
    else
        right = Set{Pair{Set{T}, Set{T}}}()
    end

    return union(node, left, right)
end