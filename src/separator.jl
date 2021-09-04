"""
    find_feg_separator_lt(skeleton, faces, root)

Finds the a valid separator for a finite element graph given its skeleton and
faces of its planar embedding represented as a set of sets of edges. The search
will initiate on a bfs tree rooted at the given vertex.
"""
function find_feg_separator_lt(skeleton::LabeledGraph{T}, faces::Set{Set{Pair{T, T}}}, root::T)::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    g = _find_finite_element_graph(skeleton, faces)

    g_star_star = copy(skeleton)
    new_vertices = Dict{Number, Set{Pair{T, T}}}()
    # Only need to add vertex to faces with more than 3 vertices   
    faces = filter(f -> length(f) > 3, faces)
    for face in faces
        new_vertex = rand(T)
        push!(new_vertices, new_vertex => face)
        lg_add_vertex!(g_star_star, new_vertex)
        face_vertices = reduce((x, y) -> union!(x, Set([y.first, y.second])), face, init=Set{Number}())
        for vertex in face_vertices
            lg_add_edge!(g_star_star, vertex, new_vertex)
        end
    end

    (c_star_star, a_star_star, b_star_star) = find_separator_lt(g_star_star, root)

    a = filter(vertex -> !(vertex in keys(new_vertices)), a_star_star)
    b = filter(vertex -> !(vertex in keys(new_vertices)), b_star_star)

    for pair in new_vertices
        if (pair.first in c_star_star)
            bad_face_vertecies =  reduce((x, y) -> union!(x, Set([y.first, y.second])), pair.second, init=Set{Number}())

            a_star_star_boundry = filter(vertex -> vertex in a_star_star, bad_face_vertecies)
            b_star_star_boundry = filter(vertex -> vertex in b_star_star, bad_face_vertecies)
            
            if length(a_star_star_boundry) < length(b_star_star_boundry)
                setdiff!(a, a_star_star_boundry)
                setdiff!(a_star_star, a_star_star_boundry)
            else
                setdiff!(b, b_star_star_boundry)
                setdiff!(b_star_star, b_star_star_boundry)
            end
        end
    end

    return (setdiff(lg_vertices(g), a, b), a, b)
end

"""
    find_feg_separator_lt_best(skeleton, faces)

Finds the a valid separator for a finite element graph given its skeleton and
faces of its planar embedding represented as a set of sets of edges. Repeats
the search from every possible BFS root and returns the most balanced
separator. 
"""
function find_feg_separator_lt_best(skeleton::LabeledGraph{T}, faces::Set{Set{Pair{T, T}}})::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    separators = map(root -> find_feg_separator_lt(skeleton, faces, root), collect(keys(skeleton.labels)))

    balances = map(separator -> min(length(separator[2]) / length(separator[3]), length(separator[3]) / length(separator[2])), separators)
    ######## See _find_feg_separator_lt_no_empty in biclique.jl
    # balances = map(separator -> begin
    #     if !isempty(separator[2]) && !isempty(separator[3])
    #         return min(length(separator[2]) / length(separator[3]), length(separator[3]) / length(separator[2]))
    #     else
    #         return -1.0
    #     end
    # end, separators)


    max_balance = max(balances...)
    for i in 1:length(balances)
        if balances[i] == max_balance return separators[i] end
    end
    return separators[1]

end

"""
    find_separator_lt(lg, root)

Finds a valid separator for a planar graph using the Lipton-Tarjan algorithm.
The search will initiate on a bfs tree rooted at the given vertex.
"""
function find_separator_lt(lg::LabeledGraph{T}, root::T)::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    levels = @pipe LightGraphs.bfs_tree(lg.graph, lg.labels[root]) |> _find_bfs_levels(_, lg.labels[root]) |> map(level -> Set(convert_vertices(lg.labels, collect(level))), _)

    # Phase I
    middle_index = 1
    for i in 2:length(levels)
        if sum(map(j -> length(levels[j]), 1:i)) > (LightGraphs.nv(lg.graph) / 2)
            middle_index = i
            break
        end
    end
    middle = levels[middle_index]

    if length(middle) < (2 * sqrt(2 * LightGraphs.nv(lg.graph)))
        above = middle_index == 1 ? Set() : @pipe map(i -> levels[i], 1:(middle_index - 1)) |> reduce(union, _)
        below = middle_index < length(levels) ? (@pipe map(i -> levels[i], (middle_index + 1):length(levels)) |> reduce(union, _)) : Set()
        a = length(above) < length(below) ? above : below;
        b = length(above) < length(below) ? below : above;
        return (middle, a, b)
    end

    # Phase II
    upper_index = -1;
    for i in middle_index:-1:1
        distance = middle_index - i
        if length(levels[i]) < 2 * (sqrt(LightGraphs.nv(lg.graph)) - distance)
            upper_index = i
            break
        end
    end
    if sum(map(i -> length(levels[i]), 1:upper_index)) > 1/3
        # upper is separator
        above = upper_index > 1 ? (@pipe map(i -> levels[i], 1:(upper_index - 1)) |> reduce(union, _)) : Set()
        below = upper_index < length(levels) ? (@pipe map(i -> levels[i], (upper_index + 1):length(levels)) |> reduce(union, _)) : Set()
        a = length(above) < length(below) ? above : below;
        b = length(above) < length(below) ? below : above;
        return (levels[upper_index], a, b)
    end

    lower_index = -1;
    for i in middle_index:length(levels)
        distance = i - middle_index
        if length(levels[i]) < 2 * (sqrt(LightGraphs.nv(lg.graph)) - distance)
            lower_index = i
            break
        end
    end
    if sum(map(i -> length(levels[i]), lower_index:length(levels))) > 1/3
        # lower is separator
        above = lower_index > 1 ? (@pipe map(i -> levels[i], 1:(lower_index - 1)) |> reduce(union, _)) : Set()
        below = lower_index < length(levels) ? (@pipe map(i -> levels[i], (lower_index + 1):length(levels)) |> reduce(union, _)) : Set()
        a = length(above) < length(below) ? above : below;
        b = length(above) < length(below) ? below : above;
        return (levels[lower_index], a, b)
    end

    # Phase III
    return find_separator_fcs(lg, root)

end

"""
    find_separator_fcs(lg, root)

Finds a valid separator for a planar graph using finite cycles as separators.
The search will initiate on a bfs tree rooted at the given vertex.
"""
function find_separator_fcs(lg::LabeledGraph{T}, root::T)::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    parents = LightGraphs.bfs_parents(lg.graph, lg.labels[root])

    tree_edges = map(tup -> LightGraphs.Edge(tup[1] => tup[2]), zip(parents, 1:length(parents)))
    non_tree_edge = filter(edge -> !(in(edge, tree_edges) || in(reverse(edge), tree_edges)), collect(LightGraphs.edges(lg.graph)))[1]

    fundamental_cycle = @pipe _find_fundamental_cycle(parents, non_tree_edge) |> Set(convert_vertices(lg.labels, collect(_)))
    return (fundamental_cycle, _find_partitions(lg, fundamental_cycle)...)
end

"""
    find_separator_fcs_best(lg, root)

Finds a valid separator for a planar graph using finite cycles as separators.
Repeats the search from every possible BFS root and returns the most balanced
separator.
"""
function find_separator_fcs_best(lg::LabeledGraph{T}, root::T)::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    # THIS ALGORITHM IS O(n^2) BECAUSE IT CHOOSES THE MOST BALANCED
    parents = LightGraphs.bfs_parents(lg.graph, lg.labels[root])

    tree_edges = map(tup -> LightGraphs.Edge(tup[1] => tup[2]), zip(parents, 1:length(parents)))
    non_tree_edges = filter(edge -> !(in(edge, tree_edges) || in(reverse(edge), tree_edges)), collect(LightGraphs.edges(lg.graph)))

    fundamental_cycles = @pipe map(non_tree_edge -> _find_fundamental_cycle(parents, non_tree_edge), non_tree_edges) |> map(cycle -> Set(convert_vertices(lg.labels, collect(cycle))), _)
    partitions = map(cycle -> (cycle, _find_partitions(lg, cycle)...), fundamental_cycles)
    # TODO: What happens when sets are empty? Can get Inf
    ######## See _find_feg_separator_lt_no_empty in biclique.jl
    balances = map(partition -> min(length(partition[2]) / length(partition[3]), length(partition[3]) / length(partition[2])), partitions)

    max_balance = max(balances...)
    for i in 1:length(balances)
        if balances[i] == max_balance return partitions[i] end
    end
    return partitions[1]
end

"""
    _find_partitions(lg, separator)

Given a graph and separator find the sets of vertices separated.
"""
function _find_partitions(lg::LabeledGraph{T}, separator::Set{T})::Tuple{Set{T}, Set{T}} where {T}
    lg_without_separator = copy(lg)

    for vertex in separator
        lg_rem_vertex!(lg_without_separator, vertex)
    end

    components = @pipe LightGraphs.connected_components(lg_without_separator.graph) |>
        map(component -> convert_vertices(lg_without_separator.labels, component), _) |>
        sort(_, by=length, rev=true)
    a = Array{T, 1}()
    b = Array{T, 1}()
    for component in components
        append!(sum(map(length, a)) < sum(map(length, b)) ? a : b, component)
    end

    return (Set(a), Set(b))
end

"""
    _find_fundamental_cycle(parents, non_tree_edge)

Given a bfs tree (in array form as provided by LightGraphs) and a non tree edge
return the fundamtal cycle as a set of vertices.
"""
function _find_fundamental_cycle(parents::Vector{Int}, non_tree_edge::LightGraphs.Edge)::Set{Int}
    left_vertex = LightGraphs.src(non_tree_edge)
    left_path = [left_vertex]
    while left_vertex != parents[left_vertex]
        left_vertex = parents[left_vertex]
        pushfirst!(left_path, left_vertex)
    end

    right_vertex = LightGraphs.dst(non_tree_edge)
    right_path = [right_vertex]
    while right_vertex != parents[right_vertex]
        right_vertex = parents[right_vertex]
        pushfirst!(right_path, right_vertex)
    end

    for i in 1:max(length(left_path), length(right_path))
        if left_path[i] != right_path[i]
            return Set(union(left_path[i:end], right_path[i:end], left_path[i - 1]))
        end
    end

    return Set(union(left_path, right_path))
end

"""
    _find_bfs_levels(tree, root)

Run breadth first search on a graph and return its BFS tree as an array of sets
representing levels of the tree.
"""
function _find_bfs_levels(tree::LightGraphs.AbstractGraph, root::Int)::Array{Set{Int}, 1}
    levels = Dict{Int, Set{Int}}()

    levels[1] = Set([root])
    i = 2
    while true
        levels[i] = Set(reduce(union, map(vertex -> LightGraphs.outneighbors(tree, vertex), collect(levels[i - 1]))))

        if isempty(levels[i]) break end

        i += 1
    end

    res = Array{Set{Int}}(undef, i - 1)
    for i in 1:length(res)
        res[i] = levels[i]
    end

    return res
end

"""
    _is_valid_separator(lg, separator, a, b)

Checks if a separator is valid on graph lg and separates A and B
"""
function _is_valid_separator(lg::LabeledGraph{T}, separator::Set{T}, a::Set{T}, b::Set{T})::Bool where T
    is_separating = true
    for source in a
        break_early = false
        for destination in b
            is_path = ClutteredEnvPathOpt.LightGraphs.has_path(lg.graph, source, destination, exclude_vertices=collect(separator))
            is_valid = is_separating && !is_path
            if !is_valid
                break_early = true
                break
            end
        end

        if break_early break end
    end

   is_covering =  all(map(vertex -> in(vertex, union(separator, a, b)), collect(keys(lg.labels))))

   is_disjoint =
    all(map(vertex -> !in(vertex, union(a, b)), collect(separator))) &&
    all(map(vertex -> !in(vertex, union(b, separator)), collect(a))) &&
    all(map(vertex -> !in(vertex, union(a, separator)), collect(b)))

    return is_separating && is_covering && is_disjoint
end

"""
    _find_finite_element_graph(skeleton, faces)

Creates a finite element graph from a skeleton and set of sets of edges
representing faces.
"""
function _find_finite_element_graph(skeleton::LabeledGraph{T}, faces::Set{Set{Pair{T, T}}})::LabeledGraph{T} where {T}
    g = copy(skeleton)
    
    for face in faces
        face_vertices = reduce((x, y) -> union!(x, Set([y.first, y.second])), face, init=Set{Number}())
        for edge in face
            for dest in face_vertices
                if (edge.first != dest)
                    lg_add_edge!(g, edge.first, dest)
                end
            end
        end
    end

    return g
end

"""
    pp_expell(lg, separator, a, b)

A postprocessing algorithm to shrink the size of a separator by expelling
vertices in the separator not connected to both A and B.
"""
function pp_expell(lg::LabeledGraph{T}, separator::Set{T}, a::Set{T}, b::Set{T})::Tuple{Set{T}, Set{T}, Set{T}} where {T}
    pp_separator = copy(separator)
    pp_a = copy(a)
    pp_b = copy(b)
    
    for separator_vertex in separator
        source = lg.labels[separator_vertex]

        is_connected_a = any(
            map(
                destination -> LightGraphs.has_path(lg.graph, source, destination, exclude_vertices=collect(convert_vertices_rev(lg.labels, setdiff(pp_separator, Set([separator_vertex]))))),
                collect(convert_vertices_rev(lg.labels, pp_a))
            )
        )

        is_connected_b = any(
            map(
                destination -> LightGraphs.has_path(lg.graph, source, destination, exclude_vertices=collect(convert_vertices_rev(lg.labels, setdiff(pp_separator, Set([separator_vertex]))))),
                collect(convert_vertices_rev(lg.labels, pp_b))
            )
        )

        if !is_connected_a && !is_connected_b
            push!(length(pp_a) < length(pp_b) ? pp_a : pp_b, separator_vertex)
            delete!(pp_separator, separator_vertex)
        elseif is_connected_a && !is_connected_b
            push!(pp_a, separator_vertex)
            delete!(pp_separator, separator_vertex)
        elseif !is_connected_a && is_connected_b
            push!(pp_b, separator_vertex)
            delete!(pp_separator, separator_vertex)
        end
    end

    # Assure A and B are nonempty
    
    # triggered = false
    # if isempty(pp_a) || isempty(pp_b)
    #     println("\n#########\nA: $pp_a\nB:$pp_b\nC:$pp_separator")
    #     triggered = true
    # end

    if isempty(pp_a)
        # Move vertices from B into A
        # First, need to move vertices in B connected to every other vertex in B into C
        for b in pp_b
            if length(pp_b) <= 1
                #println("Don't want an empty B")
                break
            end
            source = lg.labels[b]
            num_edges_with_pp_b = sum(
                map(destination -> begin
                    if LightGraphs.has_edge(lg.graph, source, destination)
                        return 1
                    else
                        return 0
                    end
                end,
                    collect(convert_vertices_rev(lg.labels, setdiff(pp_b, Set(b))))
                )
            )
            if num_edges_with_pp_b == length(setdiff(pp_b, b))
                #println("Moving $b from B into C")
                push!(pp_separator, b)
                delete!(pp_b, b)
            end
        end

        # Next we move into A either
        # 1) vertices in B that aren't adjacent to any of the other vertices in B, or
        # 2) a vertex in B that is adjacent to the least amount of other vertices in B, also moving
        # those vertices it is adjacent to into C
        size_pp_b = length(pp_b)
        added = 0
        nv = length(pp_separator) + length(pp_b) + length(pp_a)
        min = nv * (nv - 1) / 2
        min_vertex = -1

        # Case 1
        for b in pp_b
            if added <= floor(size_pp_b/2) && length(pp_b) > 1 # don't want to add more than half into A and don't want B to be empty
                source = lg.labels[b]
                num_edges_with_pp_b = sum(
                    map(destination -> begin
                        if LightGraphs.has_edge(lg.graph, source, destination)
                            return 1
                        else
                            return 0
                        end
                    end,
                        collect(convert_vertices_rev(lg.labels, setdiff(pp_b, Set(b))))
                    )
                )
                if num_edges_with_pp_b == 0
                    #println("Moving $b from B into A")
                    push!(pp_a, b)
                    delete!(pp_b, b)
                    added += 1
                elseif num_edges_with_pp_b < min
                    min = num_edges_with_pp_b
                    min_vertex = b
                end
            end
        end

        # Case 2
        if added == 0
            source = lg.labels[min_vertex]
            rev = _reverse_labels(lg.labels)
            connections = map(destination -> begin
                if LightGraphs.has_edge(lg.graph, source, destination)
                    return rev[destination]
                else
                    return 0
                end
            end,
                collect(convert_vertices_rev(lg.labels, setdiff(pp_b, Set(min_vertex))))
            )

            connections = filter(c -> c != 0, connections)

            for v in connections
                push!(pp_separator, v)
                delete!(pp_b, v)
            end

            push!(pp_a, min_vertex)
            delete!(pp_b, min_vertex)
        end


    elseif isempty(pp_b)
        # Move vertices from A into B
        # First, need to move vertices in A connected to every other vertex in A into C
        for a in pp_a
            if length(pp_a) <= 1
                #println("Don't want an empty A")
                break
            end
            source = lg.labels[a]
            num_edges_with_pp_a = sum(
                map(destination -> begin
                    if LightGraphs.has_edge(lg.graph, source, destination)
                        return 1
                    else
                        return 0
                    end
                end,
                    collect(convert_vertices_rev(lg.labels, setdiff(pp_a, Set(a))))
                )
            )
            if num_edges_with_pp_a == length(setdiff(pp_a, a))
                #println("Moving $a from A into C")
                push!(pp_separator, a)
                delete!(pp_a, a)
            end
        end
        
        # Next we move into B either
        # 1) vertices in A that aren't adjacent to any of the other vertices in A, or
        # 2) a vertex in A that is adjacent to the least amount of other vertices in A, also moving
        # those vertices it is adjacent to into C
        size_pp_a = length(pp_a)
        added = 0
        nv = length(pp_separator) + length(pp_b) + length(pp_a)
        min = nv * (nv - 1) / 2
        min_vertex = -1
        
        # Case 1
        for a in pp_a
            if added <= floor(size_pp_a/2) && length(pp_a) > 1 # don't want to add more than half into B and don't want A to be empty
                source = lg.labels[a]
                num_edges_with_pp_a = sum(
                    map(destination -> begin
                        if LightGraphs.has_edge(lg.graph, source, destination)
                            return 1
                        else
                            return 0
                        end
                    end,
                        collect(convert_vertices_rev(lg.labels, setdiff(pp_a, Set(a))))
                    )
                )
                if num_edges_with_pp_a == 0
                    #println("Moving $a from A into B")
                    push!(pp_b, a)
                    delete!(pp_a, a)
                    added += 1
                elseif num_edges_with_pp_a < min
                    min = num_edges_with_pp_a
                    min_vertex = a
                end
            end
        end

        # Case 2
        if added == 0
            source = lg.labels[min_vertex]
            rev = _reverse_labels(lg.labels)
            connections = map(destination -> begin
                if LightGraphs.has_edge(lg.graph, source, destination)
                    return rev[destination]
                else
                    return 0
                end
            end,
                collect(convert_vertices_rev(lg.labels, setdiff(pp_a, Set(min_vertex))))
            )

            connections = filter(c -> c != 0, connections)

            for v in connections
                push!(pp_separator, v)
                delete!(pp_a, v)
            end

            push!(pp_b, min_vertex)
            delete!(pp_a, min_vertex)
        end

    end

    # if triggered
    #     println("\nThe new sets are \nA: $pp_a\nB:$pp_b\nC:$pp_separator\n")
    # end

    return (pp_separator, pp_a, pp_b)
end