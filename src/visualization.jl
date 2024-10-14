# Visualizing and printing

# Overload the printing for a Degopt
function Base.show(io::IO, mime::MIME"text/plain",obj::Degopt)

    (A,B,c)=get_degopt_coeffs(obj)

    println(io,"(A,B,c)::(", string(typeof(A)),",",string(typeof(B)),",",string(typeof(c)),")");
    m=size(A,1);
    M=Matrix{String}(undef,m,2*m+2)
    M.="";

    for i=1:m
        for j=1:(m+1)

            myio=IOBuffer();
            show(myio, mime, A[i,j]);
            s=String(take!(myio))
            M[i,j]=s;
            if (A[i,j]==0)
                M[i,j]="";
            end


            myio=IOBuffer();
            show(myio, mime, B[i,j]);
            s=String(take!(myio))
            M[i,j+m+1]=s;
            if (B[i,j]==0)
                M[i,j+m+1]="";
            end
        end
    end

    M_length=mapslices( x-> maximum(length.(x)), M,dims=1);
    for i=1:size(M,1)
        for j=1:m+1
            sz=M_length[j];
            print(io,lpad(M[i,j],sz));
            print(io," ");
        end
        print("| ");
        for j=1:m+1
            sz=M_length[j+m+1];
            print(io,lpad(M[i,j+m+1],sz));
            if (j<m+1)
                print(io," ");
            end
        end
        println(io);
    end
    show(io,mime,c)


end
