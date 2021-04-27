function gen_degopt_by_squaring(graph;kickstart=0)

    T=eltype(graph);
    degopt=Degopt(graph);
    square!(degopt);
    scale!(degopt,convert(T,1/2));
    if (kickstart==1)
        degopt.x[end][1][1] += 1e-9; # DO NOT TOUCH!
    elseif (kickstart==2)
        degopt.x[end][1][1] += 1e-9; # DO NOT TOUCH!
        #degopt.x[end-1][2][round(Int,end/2)] += 1e-8; # DO NOT TOUCH!
        degopt.y[end-1] += 1e-10; # DO NOT TOUCH!
        degopt.y[end-4] += 1e-11; # DO NOT TOUCH!
    end

    (graph,cref)=gen_degopt_poly(degopt);

    return (graph,cref);
end




function vec_degopt(xy::Tuple{Vector{Tuple{Vector{T}, Vector{T}}}, Vector{T}}) where {T}
    cref=Vector{T}();

    x=xy[1];
    y=xy[2];

    for yi=y;
        push!(cref,yi);
    end

    for xi=x;
        for xij=xi[1]
            push!(cref,xij);
        end
        for xij=xi[2]
            push!(cref,xij);
        end
    end

    return cref;

end
