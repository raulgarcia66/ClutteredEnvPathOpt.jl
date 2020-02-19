module ClutteredEnvPathOpt

function bfs(graph::Dict{T, Set{T}}, root::T)::Dict{T, T} where {T}
    parents = Dict{T, T}()
    visited = Set(root)
    queue = [root]

    while !isempty(queue)
        node = pop!(queue)
        for neighbor in graph[node]
            if !in(neighbor, visited)
                push!(visited, neighbor)
                pushfirst!(queue, neighbor)
                parents[neighbor] = node
            end
        end
    end

    return parents
end

bfs(graph::Dict{T, Set{T}}) where {T} = bfs(graph, collect(keys(graph))[1])

function edges(graph::Dict{T, Set{T}})::Set{Set{T}} where {T}
    edges_per_node = map(pair -> map(dest -> Set([pair.first, dest]), collect(pair.second)), collect(graph))
    return Set(reduce(union, edges_per_node))
end

function tree_edges(parents::Dict{T, T})::Set{Set{T}} where {T}
    return Set(map(pair -> Set(pair), collect(parents)))
end

function find_fundamental_cycle(parents::Dict{T, T}, non_tree_edge::Set{T})::Set{T} where {T}
    ordered_edge = collect(non_tree_edge)
    left_node = ordered_edge[1]
    right_node = ordered_edge[2]

    left_seen = [left_node]
    right_seen = [right_node]

    while true
        left_node = get(parents, left_node, nothing)
        right_node = get(parents, right_node, nothing)

        if isnothing(left_node) && isnothing(right_node)
            return Set()
        elseif left_node == right_node
            return Set(union(left_seen, right_seen, left_node))
        elseif !isnothing(left_node) && in(left_node, right_seen)
            return Set(union(left_seen, right_seen[1:(findall(x -> x == left_node, right_seen)[1])]))
        elseif !isnothing(right_node) &&  in(right_node, left_seen)
            return Set(union(right_seen, left_seen[1:(findall(x -> x == right_node, left_seen)[1])]))
        else
            if !isnothing(left_node) push!(left_seen, left_node) end
            if !isnothing(right_node) push!(right_seen, right_node) end
        end
    end
end

function find_smaller_component(graph::Dict{T, Set{T}}, separator::Set{T})::Set{T} where {T}
    non_separator_nodes = setdiff(keys(graph), separator)
    if isempty(non_separator_nodes)
        return Set()
    end
    root = collect(non_separator_nodes)[1]

    visited = Set(root)
    queue = [root]

    while !isempty(queue)
        node = pop!(queue)
        for neighbor in graph[node]
            if !in(neighbor, visited) && !in(neighbor, separator)
                push!(visited, neighbor)
                pushfirst!(queue, neighbor)
            end
        end
    end

    # return the smaller connected component
    if length(visited) < length(graph) - (length(separator) + length(visited))
        return visited
    else
        return setdiff(keys(graph), separator, visited)
    end
end

end # module
