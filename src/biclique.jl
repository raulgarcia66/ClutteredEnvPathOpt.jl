"""
    find_biclique_cover(skeleton, faces)

Given a finite element graph's skeleton and the faces of its planar embedding,
return the biclique cover as a set of pairs of sets. This algorithm uses a
divide and conquere approach, though it has not been optimized for parallelism
nor tail-recursion.

Note: ensure the faces are passed in as a set of clockwise-ordered vertices
"""
function find_biclique_cover(skeleton::LabeledGraph{T}, faces::Set{Vector{T}})::Set{Pair{Set{T}, Set{T}}} where {T}
    face_pairs = _find_face_pairs(faces)
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)

    if ClutteredEnvPathOpt.lg_is_complete(feg)
        return Set{Pair{Set{T}, Set{T}}}()
    end

    (C, A, B) = _find_feg_separator_lt_no_empty(skeleton, face_pairs)

    skeleton_ac, faces_ac = _find_skeleton_faces(union(A, C), skeleton, faces)
    skeleton_bc, faces_bc = _find_skeleton_faces(union(B, C), skeleton, faces)

    node = Set([A => B])
    left = find_biclique_cover(skeleton_ac, faces_ac)
    right = find_biclique_cover(skeleton_bc, faces_bc)

    @show node, left, right

    return union(node, left, right)
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

    skeleton_ac, faces_ac = ClutteredEnvPathOpt._find_skeleton_faces(union(A, C), skeleton, faces)
    skeleton_bc, faces_bc = ClutteredEnvPathOpt._find_skeleton_faces(union(B, C), skeleton, faces)

    node = A => B
    left = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton_ac, faces_ac)
    right = ClutteredEnvPathOpt.find_biclique_cover_as_tree(skeleton_bc, faces_bc)

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

    faces = @pipe map(face -> filter(vertex -> in(vertex, vertices), face), collect(old_faces)) |> Set(_)

    for face in faces
        if length(face) == 2
            ClutteredEnvPathOpt.lg_add_edge!(skeleton, face[1], face[2])
        end

        if length(face) == 3
            ClutteredEnvPathOpt.lg_add_edge!(skeleton, face[1], face[2])
            ClutteredEnvPathOpt.lg_add_edge!(skeleton, face[2], face[3])
            ClutteredEnvPathOpt.lg_add_edge!(skeleton, face[1], face[3])
        end
    end

    faces = @pipe filter(face -> length(face) > 2, collect(faces)) |> Set(_)

    return skeleton, faces
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
given root either A or B is empty.
"""
function _find_feg_separator_lt_no_empty(skeleton::LabeledGraph{T}, face_pairs::Set{Set{Pair{T, T}}})::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)
    for root in keys(skeleton.labels)
        (C, A, B) = find_feg_separator_lt(skeleton, face_pairs, root)
        (C, A, B) = pp_expell(feg, C, A, B)

        if !isempty(A) && !isempty(B)
            return (C, A, B)
        end
    end
    
    return (Set([]), Set([]), Set([]))
end

"""
    _is_valid_biclique_cover(lg, cover)

Tests whether or not a cover is a valid biclique cover of a griven graph.
"""
function _is_valid_biclique_cover(lg::LabeledGraph{T}, cover::Set{Pair{Set{T}, Set{T}}})::Bool where {T}
    e_bar = LightGraphs.edges(LightGraphs.complement(lg.graph))

    edges = @pipe map(pair -> _cartesian_product(pair.first, pair.second), collect(cover)) |> reduce(union!, _, init=Set{Pair{T, T}}())
    test = LightGraphs.SimpleGraph(LightGraphs.nv(lg.graph))
    rev = _reverse_labels(lg.labels)

    for edge in edges
        LightGraphs.add_edge!(test, rev[edge.first], rev[edge.second])
    end

    return LightGraphs.edges(test) == e_bar
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
    _wrapper_find_feg_separator_lt_no_empty(skeleton, face_pairs, locations)
Wrapper function that finds (C, A, B) and plots the edges found between A and B.
"""
function _wrapper_find_feg_separator_lt_no_empty(skeleton::LabeledGraph{T}, face_pairs::Set{Set{Pair{T, T}}}, points::Vector{Any})::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    (C,A,B) = _find_feg_separator_lt_no_empty(skeleton, face_pairs)
    #(C,A,B) = find_feg_separator_lt_best(skeleton, face_pairs)
    # find_feg_separator_lt_best doesn't guarantee sets are nonempty; problem with bicliques if empty

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

    for edge_set in new_edges
        (x1,y1) = points[edge_set.first]
        (x2,y2) = points[edge_set.second]
        Plots.plot!([x1, x2], [y1, y2], color=:green, linestyle=:dash,title="Biclique # ")
    end

    display(Plots.plot!(legend=:false,xlims=(-.05,1.05),ylims=(-0.05,1.05)))

    return (C,A,B)
end

"""
    _wrapper_find_biclique_cover(skeleton, faces, locations)
Wrapper function that finds a biclique cover and plots the edges of each biclique.
"""
function _wrapper_find_biclique_cover(skeleton::LabeledGraph{T}, faces::Set{Vector{T}}, points::Vector{Any})::Set{Pair{Set{T}, Set{T}}} where {T}
    face_pairs = _find_face_pairs(faces)
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)

    if ClutteredEnvPathOpt.lg_is_complete(feg)
        return Set{Pair{Set{T}, Set{T}}}()
    end

    #(C, A, B) = _find_feg_separator_lt_no_empty(skeleton, face_pairs)
    (C, A, B) = _wrapper_find_feg_separator_lt_no_empty(skeleton, face_pairs, points)
    # (C, A, B), new_edges = _wrapper_find_feg_separator_lt_no_empty(skeleton, face_pairs, points, current_edges)
    # push!(edges, new_edges)
    skeleton_ac, faces_ac = _find_skeleton_faces(union(A, C), skeleton, faces)
    skeleton_bc, faces_bc = _find_skeleton_faces(union(B, C), skeleton, faces)

    node = Set([A => B])
    # left = find_biclique_cover(skeleton_ac, faces_ac)
    # right = find_biclique_cover(skeleton_bc, faces_bc)
    left = _wrapper_find_biclique_cover(skeleton_ac, faces_ac, points)
    right = _wrapper_find_biclique_cover(skeleton_bc, faces_bc, points)

    @show node, left, right

    return union(node, left, right)
end