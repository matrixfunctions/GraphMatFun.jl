load("n2000_2_2.mat");
m=50; tv=zeros(m,1);
for k=1:m
    tic;
    expm(A);
    tv(k)=toc;
end

format long
mean(tv)




addpath("/tmp/");

for k=1:m
    tic;
    exp_SID_m6(A);
    tv(k)=toc;
end

format long
mean(tv)


for k=1:m
    tic;
    exp_SID_m7(A);
    tv(k)=toc;
end

format long
mean(tv)
