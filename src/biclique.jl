"""
    find_biclique_cover(skeleton, faces)

Given a finite element graph's skeleton and the faces of its planar embedding,
return the biclique cover as a set of pairs of sets. This algorithm uses a
divide and conquere approach, though it has not been optimized for parallelism
nor tail-recursion.

Note: ensure the faces are passed in as a set of clockwise-ordered vertices
"""
function find_biclique_cover(skeleton::VertexSafeGraphs.VSafeGraph, faces::Set{Vector{Int}})::Set{Pair{Set{Int}, Set{Int}}}
    face_pairs = _find_face_pairs(faces)
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)

    if LightGraphs.ne(feg) == (LightGraphs.nv(feg) * (LightGraphs.nv(feg) - 1) / 2) # if feg is complete graph
        return Set{Pair{Set{Int}, Set{Int}}}()
    end

    (C, A, B) = _find_feg_separator_lt_no_empty(skeleton, face_pairs)

    skeleton_ac, faces_ac = _find_skeleton_faces(union(A, C), skeleton, faces)
    skeleton_bc, faces_bc = _find_skeleton_faces(union(B, C), skeleton, faces)

    return union(
        Set([A => B]),
        find_biclique_cover(skeleton_ac, faces_ac),
        find_biclique_cover(skeleton_bc, faces_bc)
    )
end

"""
    _find_skeleton_faces(vertices, old_skeleton, old_faces)

Find the skeleton and faces of a finite element subgraphgraph given a subset of
vertices, a finite element graph's skeleton and the faces of its planar
embedding. Returns a (skeleton, face vector) tuple.

Note: ensure the faces are passed in as a set of clockwise-ordered vertices
"""
function _find_skeleton_faces(vertices::Set{Int}, old_skeleton::VertexSafeGraphs.VSafeGraph, old_faces::Set{Vector{Int}})::Tuple{VertexSafeGraphs.VSafeGraph, Set{Vector{Int}}}
    skeleton = copy(old_skeleton)
    for vertex in LightGraphs.vertices(old_skeleton)
        if !in(vertex, vertices)
            LightGraphs.rem_vertex!(skeleton, vertex)
        end
    end

    faces = @pipe map(face -> filter(vertex -> in(vertex, vertices), face), collect(old_faces)) |> Set(_)

    for face in faces
        if length(face) == 2
            LightGraphs.add_edge!(skeleton, face[1], face[2])
        end

        if length(face) == 3
            LightGraphs.add_edge!(skeleton, face[1], face[2])
            LightGraphs.add_edge!(skeleton, face[2], face[3])
            LightGraphs.add_edge!(skeleton, face[1], face[3])
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
function _find_face_pairs(faces::Set{Vector{Int}})::Set{Set{LightGraphs.Edge}}
    face_pairs = Set{Set{LightGraphs.Edge}}()
    for face in faces
        pairs = Set{LightGraphs.Edge}()

        for i in 2:length(face)
            push!(pairs, LightGraphs.Edge(face[i - 1], face[i]))
        end

        push!(pairs, LightGraphs.Edge(face[end], face[1]))

        push!(face_pairs, pairs)
    end

    return face_pairs
end

"""
    _find_feg_separator_lt_no_empty(skeleton, face_pairs)

Find the separator of a finite element graph, repeating the process if for some
given root either A or B is empty.
"""
function _find_feg_separator_lt_no_empty(skeleton::VertexSafeGraphs.VSafeGraph, face_pairs::Set{Set{LightGraphs.Edge}})::Tuple{Set{Int}, Set{Int}, Set{Int}}
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)
    for root in LightGraphs.vertices(skeleton)
        (C, A, B) = find_feg_separator_lt(skeleton, face_pairs, root)
        # (C, A, B) = pp_expell(feg, C, A, B)

        if !isempty(A) && !isempty(B)
            return (C, A, B)
        end
    end
    
    return (Set{Int}(), Set{Int}(), Set{Int}())
end

"""
    _is_valid_biclique_cover(lg, cover)

Tests whether or not a cover is a valid biclique cover of a griven graph.
"""
function _is_valid_biclique_cover(vsg::VertexSafeGraphs.VSafeGraph, cover::Set{Pair{Set{Int}, Set{Int}}})::Bool
    e_bar = LightGraphs.edges(LightGraphs.complement(vsg.g))

    edges = @pipe map(pair -> _cartesian_product(pair.first, pair.second), collect(cover)) |> reduce(union!, _, init=Set{Pair{Int, Int}}())
    test = LightGraphs.SimpleGraph(LightGraphs.nv(vsg))

    for edge in edges
        LightGraphs.add_edge!(test, edge.first, edge.second)
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