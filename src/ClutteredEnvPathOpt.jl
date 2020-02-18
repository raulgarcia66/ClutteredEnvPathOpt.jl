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

end # module
