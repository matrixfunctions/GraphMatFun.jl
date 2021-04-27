n=MATSIZE;
A0=triu(tril(ones(n,n),3),-3)*1.0 +1.0*eye(n,n);
A0(4,7)=A0(4,7)+eps()*100; % Break symmetry to avoid special case code

col=COLUMN;

if (col==1)
    A0=2.5*A0/norm(A0,1);
else
    A0=5.5*A0/norm(A0,1);
end


fprintf("System: %s  Computer: %s\n",version,computer);
addpath('/tmp');
fprintf("Matrix norm: %d \n",norm(A0,1));
fprintf("Matrix size: %d x %d \n", n , n);
version -blas ; blas=ans;
fprintf("BLAS version: %s\n\n",blas);

nof_samples=10; tv=zeros(nof_samples,1);

% Convenience functions
expm_matlab=@(x) expm(x);
expmpoly_matlab=@(x) expmpol(x);

if (~(exist("expmpol")>0))
    fprintf("Function expmpol() not in PATH. Trying to download from\n");
    fprintf("    http://personales.upv.es/~jorsasma/software/expmpol.m\n");
    fprintf("and save in current directory.... ");
    try
        urlwrite("http://personales.upv.es/~jorsasma/software/expmpol.m","expmpol.m");
        fprintf("done.\n\n");
    catch
        fprintf("failed.\n")
        fprintf("Unable to download. expmpoly simulation set to identity.\n\n");
        expmpoly_matlab=@(x) x; % Dummy
    end
end



%%% START REPEATED CODE
% Code to run: NAME *******
fprintf("%20.20s: ","NAME");
A=A0;

for k=1:nof_samples
    tic;
    NAME(A);
    tv(k)=toc;
    pause(0.5);
end
fprintf("%.6f ± %.3f\n",median(tv),std(tv));
fprintf("                       timings:");
fprintf("%.3f ",tv);
fprintf("\n");
pause(3)

%%% END REPEATED CODE


fprintf("Done!\n ");