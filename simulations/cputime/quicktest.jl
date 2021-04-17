using LinearAlgebra,BenchmarkTools


function add_id1!(a,A,B)
    copy!(B,A);
    B+=a*I;
end

function diagview(A) # Creates a view of the diagonal
    #n=size(A,1); # assumes square
    #aa=view(A,:); dd=view(aa,1:(n+1):n^2);
    return view(A, diagind(A, 0))
end


function add_id2!(a,A,B)
    copy!(B,A);
    dd=diagview(B);
    dd .+= a;
end


function add_id3!(a,A,B,Imat);
    copy!(B,A);
    BLAS.axpby!(a,Imat,1.0,B);
end


A=randn(2000,2000);
Imat=Matrix{Float64}(I,2000,2000);
C1=zeros(2000,2000);
C2=zeros(2000,2000);
C3=zeros(2000,2000);

@btime add_id1!(3.0,A,C1);

@btime add_id2!(3.0,A,C2);

@btime add_id3!(3.0,A,C3,Imat);
