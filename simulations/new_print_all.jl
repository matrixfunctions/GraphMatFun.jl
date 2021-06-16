
include("new_set_lists_all.jl");

println("Printing simulation");

rho=mono.params[:rho];
println("Domain: $(rho), m=",m);
for (i,state)=enumerate(state_list)
    if (!isnothing(state))
        err=showerr(state,output=false,n=state.params[:target_n]*2)
        graphname=state.params[:graphname];
        println("$graphname $(Float64(err))");
    else
        println("x");
    end
end
