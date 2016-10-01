% Example with 2 objectives, n varialbles and n+1 constraints
% Demostration of the parameter 'eps'

n=200;

clear options;

% try also to vary these options:

  options.dual=0;      % use primal Benson (default)
% options.dual=1;      % use dual Benson
  options.variant='a'; % use classial variant
% options.variant='b'; % use Webster's variant (default)

P=[cos(n*2*pi./(1:1:n));-sin(n*2*pi./(1:1:n))];
B=[eye(n);-ones(1,n)];
b=[ones(n,1);-2*n];
Z=[2 1;1 2]'; % generating vectors of the dual cone of ordering cone

options.eps=1e-1;
[S1,Sh1,T1,PP1,PPh1,DD1]=bensolve(P,B,b,[],Z,[],options);
display(PP1);
display(PPh1);
display(DD1);

plotresult(PP1,PPh1,DD1,[],1,4);

options.eps=1e-2;
[S2,Sh2,T2,PP2,PPh2,DD2]=bensolve(P,B,b,[],Z,[],options);
display(PP2);
display(PPh2);
display(DD2);

plotresult(PP2,PPh2,DD2,[],2,4);

options.eps=1e-3;
[S3,Sh3,T3,PP3,PPh3,DD3]=bensolve(P,B,b,[],Z,[],options);
display(PP3);
display(PPh3);
display(DD3);

plotresult(PP3,PPh3,DD3,[],3,4);

options.eps=1e-4;
[S4,Sh4,T4,PP4,PPh4,DD4]=bensolve(P,B,b,[],Z,[],options);
display(PP4);
display(PPh4);
display(DD4);

plotresult(PP4,PPh4,DD4,[],4,4);







