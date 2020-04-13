using ClutteredEnvPathOpt
using Test
using Pipe

@testset "ClutteredEnvPathOpt.jl" begin
    # Write your own tests here.
end

@testset "fundamental cycle separator tests" begin
    GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
    RANGE_NODES = 4:32

    for graph_function in GENERATORS
        for i in RANGE_NODES
            lg = LabeledGraph(graph_function(i))
            (separator, a, b) = find_separator_fcs(lg, 1)

            @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
        end
    end
end

@testset "fundamental cycle best separator tests" begin
    GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
    RANGE_NODES = 4:32

    for graph_function in GENERATORS
        for i in RANGE_NODES
            lg = LabeledGraph(graph_function(i))
            (separator, a, b) = find_separator_fcs_best(lg, 1)

            @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
        end
    end
end

@testset "Lipton-Tarjan separator tests" begin
    GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
    RANGE_NODES = 4:32

    for graph_function in GENERATORS
        for i in RANGE_NODES
            lg = LabeledGraph(graph_function(i))
            (separator, a, b) = find_separator_lt(lg, 1)

            @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
        end
    end
end

@testset "fundamental cycle separator tests provided graphs" begin
    filenames = [
        "a280",
        "bier127",
        "ch130",
        "ch150",
        "d198",
        "d1291",
        "d1655",
        "d2103",
        "d493",
        "d657",
        "eil101",
        "eil51",
        "eil76",
        "fl1400",
        "fl1577",
        "fl417",
        "gil262",
        "kroA100",
        "kroA150",
        "kroA200",
        "kroB100",
        "kroB150",
        "kroB200",
        "kroC100",
        "kroE100"
    ]
    for filename in filenames
        open(("delaunay-graphs/$filename.tsp.del")) do file
            lines = readlines(file)
            graph = @pipe match(r"\d+", lines[1]) |> parse(Int, _.match) |> ClutteredEnvPathOpt.LightGraphs.SimpleGraph
    
            for line in lines[2:end]
                edge = map(rm -> parse(Int, rm.match), collect(eachmatch(r"\d+", line)))
                ClutteredEnvPathOpt.LightGraphs.add_edge!(graph, edge[1] + 1, edge[2] + 1); # n + 1 bc graph is 0-indexed
            end
            lg = LabeledGraph(graph)
            
            (separator, a, b) = find_separator_fcs(lg, 1)

            @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
        end
    end
end

@testset "fundamental cycle separator best tests provided graphs" begin
    filenames = [
        "a280",
        "bier127",
        "ch130",
        "ch150",
        "d198",
        "d1291",
        "d1655",
        "d2103",
        "d493",
        "d657",
        "eil101",
        "eil51",
        "eil76",
        "fl1400",
        "fl1577",
        "fl417",
        "gil262",
        "kroA100",
        "kroA150",
        "kroA200",
        "kroB100",
        "kroB150",
        "kroB200",
        "kroC100",
        "kroE100"
    ]
    for filename in filenames
        open(("delaunay-graphs/$filename.tsp.del")) do file
            lines = readlines(file)
            graph = @pipe match(r"\d+", lines[1]) |> parse(Int, _.match) |> ClutteredEnvPathOpt.LightGraphs.SimpleGraph
    
            for line in lines[2:end]
                edge = map(rm -> parse(Int, rm.match), collect(eachmatch(r"\d+", line)))
                ClutteredEnvPathOpt.LightGraphs.add_edge!(graph, edge[1] + 1, edge[2] + 1); # n + 1 bc graph is 0-indexed
            end
            lg = LabeledGraph(graph)
            
            (separator, a, b) = find_separator_fcs_best(lg, 1)

            @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
        end
    end
end

@testset "Lipton-Tarjan separator tests provided graphs" begin
    filenames = [
        "a280",
        "bier127",
        "ch130",
        "ch150",
        "d198",
        "d1291",
        "d1655",
        "d2103",
        "d493",
        "d657",
        "eil101",
        "eil51",
        "eil76",
        "fl1400",
        "fl1577",
        "fl417",
        "gil262",
        "kroA100",
        "kroA150",
        "kroA200",
        "kroB100",
        "kroB150",
        "kroB200",
        "kroC100",
        "kroE100"
    ]
    for filename in filenames
        open(("delaunay-graphs/$filename.tsp.del")) do file
            lines = readlines(file)
            graph = @pipe match(r"\d+", lines[1]) |> parse(Int, _.match) |> ClutteredEnvPathOpt.LightGraphs.SimpleGraph
    
            for line in lines[2:end]
                edge = map(rm -> parse(Int, rm.match), collect(eachmatch(r"\d+", line)))
                ClutteredEnvPathOpt.LightGraphs.add_edge!(graph, edge[1] + 1, edge[2] + 1); # n + 1 bc graph is 0-indexed
            end
            lg = LabeledGraph(graph)
            
            (separator, a, b) = find_separator_lt(lg, 1)

            @test ClutteredEnvPathOpt._is_valid_separator(lg, separator, a, b)
        end
    end
end

@testset "separator postprocessing tests" begin
GENERATORS = [ClutteredEnvPathOpt.LightGraphs.cycle_graph, ClutteredEnvPathOpt.LightGraphs.ladder_graph, ClutteredEnvPathOpt.LightGraphs.wheel_graph, x -> ClutteredEnvPathOpt.LightGraphs.grid([x, x])]
    RANGE_NODES = 4:32

    for graph_function in GENERATORS
        for i in RANGE_NODES
            lg = LabeledGraph(graph_function(i))
            (separator, a, b) = find_separator_fcs(lg, 1)
            (pp_separator, pp_a, pp_b) = pp_expell(lg, separator, a, b)

            @test ClutteredEnvPathOpt._is_valid_separator(lg, pp_separator, pp_a, pp_b)

            @test length(pp_separator) <= length(separator)
        end
    end
end