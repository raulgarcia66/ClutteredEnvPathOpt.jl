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

function convert_vertices_rev(labels::Dict{T, Int}, vertices) where {T}
    return map(vertex -> labels[vertex], collect(vertices))
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

function lg_is_complete(lg::LabeledGraph{T})::Bool where {T}
    num_edges = LightGraphs.ne(lg.graph)
    num_vertices = LightGraphs.nv(lg.graph)
    return num_edges == (num_vertices * (num_vertices - 1)) / 2
end

function _reverse_labels(labels::Dict{K, V})::Dict{V, K} where {K, V}
    return Dict(value => key for (key, value) in labels)
end
