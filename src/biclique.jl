# Faces must be given in their clockwise embedding
function find_biclique_cover(skeleton::LabeledGraph{T}, faces::Set{Vector{T}})::Set{Pair{Set{T}, Set{T}}} where {T}
    face_pairs = _find_face_pairs(faces)
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)

    if ClutteredEnvPathOpt.lg_is_complete(feg)
        return Set{Pair{Set{T}, Set{T}}}()
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

# Faces must be given in their clockwise embedding
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

function _find_feg_separator_lt_no_empty(skeleton::LabeledGraph{T}, face_pairs::Set{Set{Pair{T, T}}})::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    feg = ClutteredEnvPathOpt._find_finite_element_graph(skeleton, face_pairs)
    for root in keys(skeleton.labels)
        (C, A, B) = find_feg_separator_lt(skeleton, face_pairs, root)
        (C, A, B) = pp_expell(feg, C, A, B)
        # @show root, C, A, B

        if !isempty(A) && !isempty(B)
            return (C, A, B)
        end
    end

    println("OWIE")
    
    return ([], [], [])
end

function _is_valid_biclique_cover(lg::LabeledGraph{T}, cover::Set{Pair{Set{T}, Set{T}}})::Bool where {T}
    e_bar = LightGraphs.edges(LightGraphs.complement(lg.graph))

    edges = @pipe map(pair -> _cartesian_product(pair.first, pair.second), collect(cover)) |> reduce(union!, _, init=Set{Pair{T, T}}())
    test = LightGraphs.SimpleGraph(LightGraphs.nv(lg.graph))
    rev = _reverse_labels(lg.labels)

    for edge in edges
        LightGraphs.add_edge!(test, rev[edge.first], rev[edge.second])
    end

    res = LightGraphs.edges(test) == e_bar

    if !res
        println("OOP!!! ", length(LightGraphs.edges(test)) - length(e_bar))
        println(setdiff(LightGraphs.edges(test), e_bar))
        println(cover)
    end
    return res
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