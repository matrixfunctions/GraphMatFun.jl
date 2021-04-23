% remember to disable CPU threads boosting:
% echo "1" > /sys/devices/system/cpu/intel_pstate/no_turbo
% echo "0" > /sys/devices/system/cpu/intel_pstate/no_turbo
% from https://askubuntu.com/questions/619875/disabling-intel-turbo-boost-in-ubuntu

display("Note that the turbo should preferrably be disable. No-turbo parameter");
type('/sys/devices/system/cpu/intel_pstate/no_turbo')


%load("n2000_2_2.mat");
n=2000;
A0=triu(tril(ones(n,n),3),-3)*1.0 +1.0*eye(n,n);
%A0=A0-2*spdiags(ones(n,1),-1,n,n);

A0(4,7)=A0(4,7)+eps()*100;
A0=2.5*A0/norm(A0,1);
%A0=5.5*A0/norm(A0,1);
A=A0;

%norm(A)
%max(abs(eig(A)))


m=6;tv=zeros(m,1);


"expm native"
for k=1:m
    tic;
    expm(A);
    tv(k)=toc;
end


format long
mean(tv)
pause(2);


addpath("/tmp/");
"exp_m6_mono_taylor_2_7"
for k=1:m
    tic;
    exp_m6_mono_taylor_2_7(A);
    tv(k)=toc;
end

format long
mean(tv)
pause(2);

"exp_m6_SID_2_22"
for k=1:m
    tic;
    exp_m6_SID_2_22(A);
    tv(k)=toc;
end

format long
mean(tv)


"exp_m7_SIDplus_6_0"
for k=1:m
    tic;
    exp_m7_SIDplus_6_0(A);
    tv(k)=toc;
end

format long
mean(tv)



"exp_m7_mono_taylor_6_0"
for k=1:m
    tic;
    exp_m7_mono_taylor_6_0(A);
    tv(k)=toc;
end

format long
mean(tv)


"exp_native_73_jl"
for k=1:m
    tic;
    exp_native_73_jl(A);
    tv(k)=toc;
end

format long
mean(tv)
pause(2)

"exp_native_83_jl"
for k=1:m
    tic;
    exp_native_83_jl(A);
    tv(k)=toc;
end

format long
mean(tv)

"expm native"
for k=1:m
    tic;
    expm(A);
    tv(k)=toc;
end

format long
mean(tv)



"expmpol"
for k=1:m
    tic;
    expmpol(A);
    tv(k)=toc;
end

format long
mean(tv)
