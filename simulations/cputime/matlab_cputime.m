% remember to disable CPU threads boosting:
% Disable:
% $ sudo -s
% # echo "1" > /sys/devices/system/cpu/intel_pstate/no_turbo
% Enable:
% $ sudo -s
% # echo "0" > /sys/devices/system/cpu/intel_pstate/no_turbo
% from https://askubuntu.com/questions/619875/disabling-intel-turbo-boost-in-ubuntu

disp("Note that the turbo should preferrably be disabled.");
disp("  Next line = 1 => Disabled.")
disp("  Next line = 0 => Enabled.");
type('/sys/devices/system/cpu/intel_pstate/no_turbo')


% Test matrix
n=2000;
A0=triu(tril(ones(n,n),3),-3)*1.0 +1.0*eye(n,n);

A0(4,7)=A0(4,7)+eps()*100; % Break symmetry to avoid special case code
col=1;
if (col==1)
    A0=2.5*A0/norm(A0,1);
else
    A0=5.5*A0/norm(A0,1);
end

A=A0;


disp(strcat("norm(A,1)=", num2str(norm(A,1))))

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

% Assume generated code is in "/tmp"
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
