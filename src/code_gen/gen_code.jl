include("gen_code_mem.jl")
include("gen_code_snippets.jl")

include("gen_c_code.jl")
include("gen_julia_code.jl")
include("gen_matlab_code.jl")
include("gen_degopt_code.jl")

export gen_code

# Every language needs:
# comment(::Lang,s)
# slotname(::Lang,i) #
# assign_coeff(::Lang,v,i)
# function_definition(::Lang,T,funname)
# function_init(lang::Lang,T,mem,graph)
# init_mem(lang::Lang,max_nof_nodes)
# function_end(lang::Lang,graph,mem)
# execute_operation!(lang::Lang,T,graph,node,
#                    dealloc_list,   mem)

"""
    gen_code(fname,graph; priohelp=Dict{Symbol,Float64}(),
             lang=LangJulia(),funname="dummy")

Generates the code for the `graph` in the language
 specified in `lang` and writes it into the file
`fname`. The string `funname` is the function name.
Topological order of the nodes is comptued using
`get_topo_order` and `priohelp` can be used to
influence the order.

Currently supported languages: `LangC_MKL`, `LangC_OpenBLAS`,
 `LangJulia`, `LangMatlab`.

"""
function gen_code(fname,graph;
                  priohelp=Dict{Symbol,Float64}(),
                  lang=LangJulia(),
                  funname="dummy",
                  generate_main=false)

    if has_trivial_nodes(graph)
        error("Please run compress_graph!() on the graph first.")
    end
    T=eltype(eltype(typeof(graph.coeffs.vals)))

    if (fname isa String)
        file = open(abspath(fname), "w+")
    else
        # Lazy: Print out to stdout if no filename
        file=Base.stdout
    end

    (order, can_be_deallocated, max_nof_slots) =
        get_topo_order(graph; priohelp=priohelp)

    # max_nof_slots is the path width which gives
    # a bound on the number memory slots needed.

    println(file,to_string(function_definition(lang,T,funname)))

    # We do a double sweep in order to determine exactly
    # how many memory slots are needed. The first sweep
    # we carry out all operations but store only the
    # maximum number of memory slots needed. The
    # second sweep generates the code.

    mem=init_mem(lang,max_nof_slots+3)
    function_init(lang,T,mem,graph)

    # Sweep 1: Determine exactly the number of slots needed
    nof_slots=0
    for (i,node) in enumerate(order)
        (exec_code,result_variable)=execute_operation!(lang,
                                                       T,graph,node,
                                                       can_be_deallocated[i],
                                                       mem)

        # How many slots needed to reach this point
        if (!isnothing(findlast(mem.slots .!= :Free)))
            nof_slots=max(nof_slots,findlast(mem.slots .!= :Free))
        end
    end

    # Sweep 2:
    mem=init_mem(lang,nof_slots)
    println(file,to_string(function_init(lang,T,mem,graph)))

    println(file,comment(lang,
                         "Computation order: "*join(string.(order)," ")))
    for (i,node) in enumerate(order)
        (exec_code,result_variable)=execute_operation!(lang,
                                                       T,graph,node,
                                                       can_be_deallocated[i],
                                                       mem)
        println(file,to_string(exec_code))
    end
    println(file,to_string(function_end(lang,graph,mem)))

    # Generate main function, if necessary.
    if generate_main && typeof(lang) <: LangC
        exec_code=gen_main(lang,T,fname,funname)
        println(file,to_string(exec_code))
    end

    if (fname isa String)
        close(file)
    end

end
