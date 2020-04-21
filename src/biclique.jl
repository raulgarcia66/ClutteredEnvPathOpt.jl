function find_biclique_cover(skeleton::LabeledGraph{T}, faces::Set{Set{Pair{T, T}}})::Set{Pair{Set{T}, Set{T}}} where {T}
    res = Set()

    biclique_tree = _find_biclique_cover_helper(skeleton, faces)
    flat = flatten(biclique_tree)

    for i in 2:2:length(flat)
        if !isnothing(flat[i]) push!(res, Pair(flat[i], flat[i + 1])) end
    end

    return res
end

function _find_biclique_cover_helper(skeleton::LabeledGraph{T}, faces::Set{Set{Pair{T, T}}})::BinaryTree{Set{T}} where {T}
    if lg_is_complete(skeleton) || LightGraphs.ne(skeleton.graph) == 0
        return BinaryTree(lg_vertices(skeleton), nothing, nothing)
    end

    # compute separator
    (C, A, B) = find_feg_separator_lt(skeleton, faces, collect(keys(skeleton.labels))[1])
    println(A, B, C)

    if isempty(A)
        A = C
    end

    if isempty(B)
        B = C
    end

    skeleton_ac = copy(skeleton)
    for vertex in B
        ClutteredEnvPathOpt.lg_rem_vertex!(skeleton_ac, vertex)
    end
    faces_ac = filter(face -> issubset(reduce((x, y) -> union!(x, Set([y.first, y.second])), face, init=Set{Number}()), union(A, C)), copy(faces))

    skeleton_bc = copy(skeleton)
    for vertex in A
        ClutteredEnvPathOpt.lg_rem_vertex!(skeleton_bc, vertex)
    end
    faces_bc = filter(face -> issubset(reduce((x, y) -> union!(x, Set([y.first, y.second])), face, init=Set{Number}()), union(B, C)), copy(faces))

    
    return BinaryTree(lg_vertices(skeleton), _find_biclique_cover_helper(skeleton_ac, faces_ac), _find_biclique_cover_helper(skeleton_bc, faces_bc))
end

function _is_valid_biclique_cover(lg::LabeledGraph{T}, cover::Set{Pair{Set{T}, Set{T}}})::Bool where {T}
    println("Cover:")
    for pair in cover
        println(pair)
    end
    e_bar = LightGraphs.edges(LightGraphs.complement(lg.graph))

    edges = @pipe map(pair -> _cartesian_product(pair.first, pair.second), collect(cover)) |> reduce(union!, _, init=Set{Pair{T, T}}())
    test = LightGraphs.SimpleGraph(LightGraphs.nv(lg.graph))
    rev = _reverse_labels(lg.labels)

    for edge in edges
        LightGraphs.add_edge!(test, rev[edge.first], rev[edge.second])
    end

    println("Test edges: ", collect(LightGraphs.edges(test)))
    println("E bar: ", collect(e_bar))
    return LightGraphs.edges(test) == e_bar
end

function _cartesian_product(A::Set{V}, B::Set{W})::Set{Pair{V, W}} where {V, W}
    res = Set{Pair{V, W}}()

    for a in A
        for b in B
            push!(res, Pair(a, b))
        end
    end

    return res
end