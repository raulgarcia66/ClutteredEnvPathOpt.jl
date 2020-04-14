struct LabeledGraph{T}
    graph::LightGraphs.AbstractGraph
    labels::Dict{T, Int}
end

function LabeledGraph(g::LightGraphs.AbstractGraph)
    return LabeledGraph(g, Dict(i => i for i = LightGraphs.vertices(g)))
end

Base.copy(lg::LabeledGraph) = LabeledGraph(copy(lg.graph), copy(lg.labels))

function lg_vertices(lg::LabeledGraph{T})::Set{T} where {T}
    return Set(keys(lg.labels))
end

function convert_vertices(labels::Dict{T, Int}, vertices) where {T}
    rev = _reverse_labels(labels)
    return map(vertex -> rev[vertex], vertices)
end

function lg_rem_vertex!(lg::LabeledGraph{T}, v::T)::Bool where {T}
    rev = _reverse_labels(lg.labels)

    n = lg.labels[v]
    last = rev[LightGraphs.nv(lg.graph)]
    lg.labels[last] = n
    delete!(lg.labels, v)
    return LightGraphs.rem_vertex!(lg.graph, n)
end

function lg_add_vertex!(lg::LabeledGraph{T}, v::T)::Bool where {T}
    lg.labels[v] = LightGraphs.nv(lg.graph) + 1
    return LightGraphs.add_vertex!(lg.graph)
end

function lg_add_edge!(lg::LabeledGraph{T}, x::T, y::T)::Bool where {T}
    return LightGraphs.add_edge!(lg.graph, lg.labels[x], lg.labels[y])
end

function _reverse_labels(labels::Dict{K, V})::Dict{V, K} where {K, V}
    return Dict(value => key for (key, value) in labels)
end

struct BinaryTree{T}
    value::T
    left::Union{BinaryTree{T}, Nothing}
    right::Union{BinaryTree{T}, Nothing}
end

Base.length(tree::BinaryTree{T}) where {T} = begin
    left_length = !isnothing(tree.left) ? length(tree.left) : 0
    right_length = !isnothing(tree.right) ? length(tree.right) : 0

    return 1 + left_length + right_length
end

function height(tree::BinaryTree{T}) where {T}
    left_height = !isnothing(tree.left) ? length(tree.left) : 0
    right_height = !isnothing(tree.right) ? length(tree.right) : 0

    return 1 + max(left_height, right_height)
end

function flatten(tree::BinaryTree{T}, array::Array{Union{T, Nothing}, 1}, index::Int) where {T}
    array[index] = tree.value

    if !isnothing(tree.left) flatten(tree.left, array, index * 2) end
    if !isnothing(tree.right) flatten(tree.right, array, index * 2 + 1) end
end

function flatten(tree::BinaryTree{T})::Array{Union{T, Nothing}, 1} where {T}
    array = Array{Union{T, Nothing}, 1}(undef, height(tree) ^ 2 - 1)
    flatten(tree, array, 1)
    return array
end
