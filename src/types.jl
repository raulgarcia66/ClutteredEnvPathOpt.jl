# Overriding some LightGraphs methods to be compatible with VertexSafeGraphs.VSafeGraph

"""
    has_path(g::VSafeGraph, u, v; exclude_vertices=Vector())
Return `true` if there is a path from `u` to `v` in `g` (while avoiding vertices in
`exclude_vertices`) or `u == v`. Return false if there is no such path or if `u` or `v`
is in `excluded_vertices`. 
"""
function LightGraphs.has_path(vsg::VertexSafeGraphs.VSafeGraph{T}, u::Integer, v::Integer; 
        exclude_vertices::AbstractVector = Vector{T}()) where T
    seen = zeros(Bool, LightGraphs.nv(vsg.g))
    for ve in exclude_vertices # mark excluded vertices as seen
        seen[ve] = true
    end
    (seen[u] || seen[v]) && return false
    u == v && return true # cannot be separated
    next = Vector{T}()
    push!(next, u)
    seen[u] = true
    while !isempty(next)
        src = popfirst!(next) # get new element from queue
        for vertex in LightGraphs.outneighbors(vsg.g, src)
            vertex == v && return true
            if !seen[vertex]
                push!(next, vertex) # push onto queue
                seen[vertex] = true
            end
        end
    end
    return false

    
end

"""
    connected_components(vsg)
Return the [connected components](https://en.wikipedia.org/wiki/Connectivity_(graph_theory))
of an undirected graph `g` as a vector of components, with each element a vector of vertices
belonging to the component. When running on a VSafeGraph, LightGraphs implementation will
return deleted vertices as components of only one vertex; we must filter them out here.
"""
function LightGraphs.connected_components(vsg::VertexSafeGraphs.VSafeGraph)
    return @pipe LightGraphs.connected_components(vsg.g) |> filter(component -> !(first(component) in vsg.deleted_vertices), _)
end

function LightGraphs._bfs_parents(vsg::VertexSafeGraphs.VSafeGraph{T}, source, neighborfn::Function) where T
    n = LightGraphs.nv(vsg.g)
    visited = falses(n)
    parents = zeros(T, n)
    cur_level = Vector{T}()
    sizehint!(cur_level, n)
    next_level = Vector{T}()
    sizehint!(next_level, n)
    @inbounds for s in source
        visited[s] = true
        push!(cur_level, s)
        parents[s] = s
    end
    while !isempty(cur_level)
        @inbounds for v in cur_level
            @inbounds @simd for i in  neighborfn(vsg, v)
                if !visited[i]
                    push!(next_level, i)
                    parents[i] = v
                    visited[i] = true
                end
            end
        end
        empty!(cur_level)
        cur_level, next_level = next_level, cur_level
        sort!(cur_level)
    end
    return parents
end