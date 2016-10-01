% Example with 2 objectives, n varialbles and n+1 constraints with
% different settings

n=100;

clear options;
options.info=2;        % show more information
options.variant='b';   % use Webster's variant of Benson's algorithm
options.vert_enum='A'; % use vertex enumeration variant A
 
P=[cos(n*2*pi./(1:1:n));-sin(n*2*pi./(1:1:n))];
B=[eye(n);-ones(1,n)];
b=[ones(n,1);-2*n];
Z=[2 1;1 2]';         % generating vectors of dual cone of ordering cone


options.lp_solver=1;  % use cdd LP solver (criss-cross method)
options.dual=0;       % use primal Benson
tic
[S1,Sh1,T1,PP1,PPh1,DD1]=bensolve(P,B,b,[],Z,[],options);
toc

options.lp_solver=1;  % use cdd LP solver (criss-cross method)
options.dual=1;       % use dual Benson
tic
[S2,Sh2,T2,PP2,PPh2,DD2]=bensolve(P,B,b,[],Z,[],options);
toc

options.lp_solver=3;  % use glpk LP solver (revised simplex)
options.dual=0;       % use primal Benson
tic
[S3,Sh3,T3,PP3,PPh3,DD3]=bensolve(P,B,b,[],Z,[],options);
toc

options.lp_solver=3;  % use glpk LP solver (revised simplex)
options.dual=1;       % use dual Benson
tic
[S4,Sh4,T4,PP4,PPh4,DD4]=bensolve(P,B,b,[],Z,[],options);
toc

display(PP1);
display(PP2);
display(PP3);
display(PP4);






