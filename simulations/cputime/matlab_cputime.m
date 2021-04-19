% remember to disable CPU threads boosting:
% echo "1" > /sys/devices/system/cpu/intel_pstate/no_turbo
% echo "0" > /sys/devices/system/cpu/intel_pstate/no_turbo
% from https://askubuntu.com/questions/619875/disabling-intel-turbo-boost-in-ubuntu

%load("n2000_2_2.mat");
n=2000;
A0=triu(tril(ones(n,n),3),-3)*1.0 +1.0*eye(n,n);
%A0=A0-2*spdiags(ones(n,1),-1,n,n);

A0(4,7)=A0(4,7)+eps()*100;
A0=2.5*A0/norm(A0,1);
A=A0;

%norm(A)
%max(abs(eig(A)))


m=6;tv=zeros(m,1);


"expm native"
for k=1:m
    k
    tic;
    expm(A);
    tv(k)=toc;
end

format long
mean(tv)



addpath("/tmp/");
"exp_m6_mono_taylor_2_7"
for k=1:m
    tic;
    exp_m6_mono_taylor_2_7(A);
    tv(k)=toc;
end

format long
mean(tv)

"exp_m6_SID_2_22"
for k=1:m
    tic;
    exp_m6_SID_2_22(A);
    tv(k)=toc;
end

format long
mean(tv)


"exp_m7_SIDplus_3_59"
for k=1:m
    tic;
    exp_m7_SIDplus_3_59(A);
    tv(k)=toc;
end

format long
mean(tv)

"expm native"
for k=1:m
    k
    tic;
    expm(A);
    tv(k)=toc;
end

format long
mean(tv)
