# Read and write cgs files
using Dates
export export_compgraph, import_compgraph;

"""
    function export_compgraph(
            graph,
        fname;
        main_output = nothing,
        order = get_topo_order(graph)[1],
        fun = "",
        dom = "",
        err = "",
        genby = "",
        descr = "",
        user = splitdir(homedir())[end],
    )

Exports the graph to the file `fname` using the computation graph format (cgr).
These kwargs are strings or values that are stored as comments in the file

  - `descr` of the graph
  - `fun` function it approximates
  - `dom` domain it approximates
  - `err` error in this domain
  - `genby` script that generated this file
  - `user` a username that created the file

These kwarg influence how the graph is stored:

  - `order` specifies an order the graph nodes are stored
  - `main_output` is a `Symbol` specifying what is the output.

CGR format:

Every line corresponds to a node / operation. The syntax is matlab compatible,
although [`gen_code`](@ref) produces faster matlab code.

See also [`export_compgraph`](@ref).
"""
function export_compgraph(
    graph,
    fname;
    main_output = nothing,
    order = get_topo_order(graph)[1],
    fun = "",
    dom = "",
    err = "",
    genby = "",
    descr = "",
    user = splitdir(homedir())[end],
)
    # USER extraction for windows and linux:
    # https://stackoverflow.com/questions/27809812/get-current-username-in-julia-linux
    fname = abspath(fname)
    file = open(fname, "w+")

    nof_written_outputs = 0
    T = eltype(graph)

    println(file, "%# Representation of a computation graph")
    if !isempty(fun)
        println(file, "%# Function: ", fun)
    end
    if !isempty(descr)
        println(file, "%# ", descr)
    end
    if !isempty(dom)
        println(file, "%# Domain: ", dom)
    end
    if !isempty(err)
        println(file, "%# Error: ", err)
    end
    if !isempty(genby)
        println(file, "%# Generated by command: ", genby)
    end

    println(file, "%# Created: ", now(), " by user ", user)
    println(file, "")
    println(file, "graph_coeff_type=\"", string(T), "\";")
    println(file, "")
    for node in order
        op = graph.operations[node]
        if op == :mult
            println(
                file,
                String(node),
                "=",
                String(graph.parents[node][1]),
                "*",
                String(graph.parents[node][2]),
                ";",
            )
        elseif op == :ldiv
            println(
                file,
                String(node),
                "=",
                String(graph.parents[node][1]),
                "\\",
                String(graph.parents[node][2]),
                ";",
            )
        elseif op == :lincomb
            for (i,c) = enumerate(graph.coeffs[node])
                print(file, "coeff$i=", real(c))
                if T <: Complex
                    t = imag(c)
                    if t >= 0
                        print(file, " + ", t)
                    else
                        print(file, " - ", abs(t))
                    end
                    print(file, "i")
                end
                println(file, ";")
            end
            terms=map(x-> "coeff$(x[1])*$(x[2])",
                      enumerate(String.(graph.parents[node])))
            println(
                file,
                String(node),
                "=",
                join(terms,"+"),
                ";",
            )
        end

        if any(graph.outputs .== node) && !(node == main_output)
            println(file, "output$nof_written_outputs", "=", String(node))
            nof_written_outputs += 1
        end
    end

    if !isnothing(main_output)
        println(file, "output$nof_written_outputs", "=", String(main_output))
    end
    return close(file)
end

"""
    graph=import_compgraph(fname)

Reads a graph stored in the computation graph format (cgr) in the file `fname`.

See also [`import_compgraph`](@ref).
"""
function import_compgraph(fname)
    fname = abspath(fname)
    file = open(fname, "r")
    line = readline(file)
    while !eof(file) && !startswith(line, "graph_coeff_type=\"")
        line = readline(file)
    end
    num = split(line, r"^.*=")[2] #remove everything upto and including "="
    num = num[2:end-2] #remove trailing semicolon and quotation marks
    T = eval(Meta.parse(num))
    graph = Compgraph(T)

    while !eof(file)
        line = readline(file)

        if startswith(line, "coeff1=") #Start of three lines defining :lincomb
            α = Vector{T}(undef, 2)
            for i = 1:2
                # Remove everything upto "=" inclusive.
                num = replace(line, r".*=" => "")
                α[i] = parse(T, num[1:end-1])
                line = readline(file)
            end
            parts = split(line, "=") # Split string on "="
            target = Symbol(parts[1])
            subparts = split(parts[2], "+")
            if endswith(subparts[2], ";")
                subparts[2] = subparts[2][1:end-1]
            end
            p = Vector{Symbol}(undef, 2)
            for i = 1:2
                # Remove everything upto "*" inclusive
                p[i] = Symbol(replace(subparts[i], r".*\*" => ""))
            end
            add_lincomb!(graph, target, α[1], p[1], α[2], p[2])
        elseif occursin(r"\*", line) # Looks for "*". Found :mult
            parts = split(line, "=") # Split string on "="
            target = Symbol(parts[1])
            subparts = split(parts[2], "*")
            if endswith(subparts[2], ";")
                subparts[2] = subparts[2][1:end-1]
            end
            p = Vector{Symbol}(undef, 2)
            for i = 1:2
                p[i] = Symbol(subparts[i])
            end
            add_mult!(graph, target, p[1], p[2])

        elseif occursin(r"\\\\", line) #Looks for "\". Found :ldiv
            parts = split(line, "=") # Split string on "="
            target = Symbol(parts[1])
            subparts = split(parts[2], "\\")
            if endswith(subparts[2], ";")
                subparts[2] = subparts[2][1:end-1]
            end
            p = Vector{Symbol}(undef, 2)
            for i = 1:2
                p[i] = Symbol(subparts[i])
            end
            add_ldiv!(graph, target, p[1], p[2])

        elseif startswith(line, r"output[0-9][0-9]*=") #Define output
            parts = split(line, "=") # Split string on "="
            out = Symbol(parts[2])
            add_output!(graph, out)
        end
    end

    close(file)
    return graph
end
