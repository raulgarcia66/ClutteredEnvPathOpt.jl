function find_biclique_cover(lg:LabeledGraph{T})::Set{Pair{Set{T}, Set{T}}} where {T}
    res = Set()

    biclique_tree = _find_biclique_cover_helper(lg)
    flat = flatten(biclique_tree)

    for i in 2:2:length(flat)
        push!(res, Pair(flat[i], flat[i + 1]))
    end

    return res
end

function _find_biclique_cover_helper(lg:LabeledGraph{T})::BinaryTree{T} where {T}

end