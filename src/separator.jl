using LightGraphs

function fundamental_cycle_separator(graph::SimpleGraph, root::Int)::Set{Int}
    parents = bfs_parents(graph, root)

    tree_edges = map(tup -> Edge(tup[1] => tup[2]), zip(parents, 1:length(parents)))
    non_tree_edges = filter(edge -> !(in(edge, tree_edges) || in(reverse(edge), tree_edges)), collect(edges(graph)))

    fundamental_cycles = map(non_tree_edge -> _find_fundamental_cycle(parents, root, non_tree_edge), non_tree_edges)

    balances = map(separator -> _find_balance(graph, Set(separator)), fundamental_cycles)
    max_balance = max(balances...)
    for i in 1:length(balances)
        if balances[i] == max_balance return fundamental_cycles[i] end
    end
    return fundamental_cycles[1]
end
fundamental_cycle_separator(graph::SimpleGraph) = fundamental_cycle_separator(graph, 1)

function _find_fundamental_cycle(parents::Array{Int, 1}, root::Int64, non_tree_edge::Edge)::Set{Int}
    left_vertex = src(non_tree_edge)
    left_path = [left_vertex]
    while left_vertex != parents[left_vertex]
        left_vertex = parents[left_vertex]
        push!(left_path, left_vertex)
    end

    right_vertex = dst(non_tree_edge)
    right_path = [right_vertex]
    while right_vertex != parents[right_vertex]
        right_vertex = parents[right_vertex]
        push!(right_path, right_vertex)
    end

    for i in 1:max(length(left_path), length(right_path))
        if left_path[i] != right_path[i]
            return Set(union(left_path[i:end], right_path[i:end], left_path[i]))
        end
    end

    return Set(union(left_path, right_path))
end

function _find_balance(graph::SimpleGraph, separator::Set{Int})::Number
    separated_edges = filter(edge -> !(in(src(edge), separator) || in(dst(edge), separator)), collect(edges(graph)))
    separated = SimpleGraphFromIterator(separated_edges)

    # filtration must happen because LightGraphs leaves in the separator nodes as orphans
    components = filter(component -> !(length(component) == 1 && in(component[1], separator)), connected_components(separated))

    if length(components) < 2
        return 0
    else
        return min(length(components[1]) / length(components[2]), length(components[2]) / length(components[1]))
    end
end