#import GraphMatFun.function_definition
#import GraphMatFun.comment
#import GraphMatFun.to_string
#import GraphMatFun.init_code
#import GraphMatFun.push_code!
#import GraphMatFun.gen_code
export LangDegoptJulia

struct LangDegoptJulia
    lang::LangJulia
end
"""
    LangDegoptJulia(lang)

Code generation specifically optimized for `Degopt` objects. Language
specific parameters are given in the underlying `lang` object.

"""
function LangDegoptJulia()
    return LangDegoptJulia(LangJulia());
end

comment(x::LangDegoptJulia,s)=comment(x.lang,s);

# Code generation.
function_definition(lang::LangDegoptJulia,graph,T,funname)=
    function_definition(lang.lang,graph,T,funname);


function function_end(lang::LangDegoptJulia,graph,mem)
    code=init_code(lang)
    push_code!(code,"return Y") # Output
    push_code!(code,"end");
    return code;
end

function fix_I_op(str::String)
    return replace(str,"I .+ " => "I + ");
end


function _gen_code(fname,graph,
                  lang::LangDegoptJulia,
                  priohelp,
                  funname)
    degopt=Degopt(graph);
    T=eltype(eltype(typeof(graph.coeffs.vals)))

    if (fname isa String)
        file = open(abspath(fname), "w+")
    else
        # Lazy: Print out to stdout if no filename
        file=Base.stdout
    end

    println(file,to_string(function_definition(lang.lang,graph,T,funname)))


    code=init_code(lang)
    #
    x_varnames=[:I,:A];
    for j=1:10
        push!(x_varnames,Symbol("X$j"));
    end

    x=degopt.x;


    coefflist=[];
    for i=1:size(x,1)
        for p=1:2
            xip_terms=Vector();
            for j=1:size(degopt.x[i][p],1)
                c=x[i][p][j];
                if (c!=0)
                    if (c==1)
                        push!(xip_terms,"$(x_varnames[j])");
                    else
                        push!(coefflist,c);
                        q=size(coefflist,1);
                        push!(xip_terms,"coeff[$q].*$(x_varnames[j])");
                    end
                end

            end
            xip=fix_I_op(join(xip_terms," .+ "));
            push_code!(code,"$(x_varnames[i+2])_$p=$xip");
        end
        push_code!(code,"$(x_varnames[i+2])=$(x_varnames[i+2])_1 * $(x_varnames[i+2])_2");
    end
    y=degopt.y;
    y_terms=Vector();
    for i=1:size(y,1)
        c=y[i]
        if (c!=0)
            if (c==1)
                push!(y_terms,"$(x_varnames[i])");
            else
                push!(coefflist,c);
                q=size(coefflist,1);
                push!(y_terms,"coeff[$q].*$(x_varnames[i])");
            end
        end

    end
    push_code!(code,"Y="*fix_I_op(join(y_terms," .+ ")));


    println(file,"coeff=T["*join(coefflist,";")*"]");
    println(file,to_string(code))


    code=function_end(lang,graph,Nothing)

    println(file,to_string(code))

    if (fname isa String)
        close(file)
    end


end
#
#
#(g,_)=gen_sid_exp(7);
#lang=LangDegoptJulia();
#gen_code("/tmp/asd.jl",g,lang);
#
#include("/tmp/asd.jl");
#A=randn(3,3);
#@show norm(dummy(A)-eval_graph(g,A))
#
#
#gen_code("/tmp/asd2.jl",g,lang=LangJulia(),funname="dummy2");
#
#include("/tmp/asd2.jl");
#@show norm(dummy2(deepcopy(A))-eval_graph(g,A))
#
#
#n=2000;
#A=randn(n,n);
#@btime dummy(A);
#@btime dummy2(A);
#
