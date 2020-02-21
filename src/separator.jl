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

function pp_expell(graph::SimpleGraph, separator::Set{Int})::Set{Int}
    res = copy(separator)
    component_a = _find_component(graph, separator)
    component_b = Set(setdiff(vertices(graph), component_a, separator))

    # if |a| == 0 then first fill it with serparator vertices not adjacent to b
    if isempty(component_a)
        for vertex in res
            if isempty(intersect(neighbors(graph, vertex), component_b))
                delete!(res, vertex)
                push!(component_a, vertex)
                break
            end
        end

        # if a started empty and nothing could be initially added expullsion cannot occur
        if isempty(component_a) return separator end
    end

    if isempty(component_b)
        for vertex in res
            if isempty(intersect(neighbors(graph, vertex), component_a))
                delete!(res, vertex)
                push!(component_b, vertex)
                break
            end
        end

        # if b started empty and nothing could be initially added expullsion cannot occur
        if isempty(component_b) return separator end 
    end

    changes_made = true
    while changes_made
        changes_made = false
        for vertex in res
            adjacent_to_a = !isempty(intersect(neighbors(graph, vertex), component_a))
            adjacent_to_b = !isempty(intersect(neighbors(graph, vertex), component_b))
            
            if (!adjacent_to_a && adjacent_to_b)
                changes_made = true
                delete!(res, vertex)
                push!(component_b, vertex)
            elseif (adjacent_to_a && !adjacent_to_b)
                changes_made = true
                delete!(res, vertex)
                push!(component_a, vertex)
            end
        end
    end

    return res
end

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

function _find_component(graph::SimpleGraph, separator::Set{Int})::Set{Int}
    non_seporator_vertices = Set(setdiff(vertices(graph), separator))
    if isempty(non_seporator_vertices)
        return Set{Int}()
    end

    root = collect(non_seporator_vertices)[1]
    component = Set(filter(vertex -> has_path(graph, root, vertex, exclude_vertices=collect(separator)), non_seporator_vertices))

    # return component
    # Return the smaller component
    if 2 * length(component) < nv(graph) - length(separator)
        return component
    end
    return setdiff(non_seporator_vertices, component)
end