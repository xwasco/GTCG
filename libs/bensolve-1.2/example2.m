% Example with 2 objectives

B=[[3 1];[1 2];[1 1]];
b=[0; 0; 1];
P=[1 0 ; 0 1];

% generating vectors of ordering cone C
Y=[-1 3; 3/2 -1];

% geometric duality parameter
c=[0;1];

clear options;
options.dual=1;        % use dual Benson algorithm
options.lp_solver=3;   % use glpk (revised simplex)
options.eps=1e-4;      % set eps used by Algorithms 1 and 2
options.variant='b';   % use Webster's variant of Algorithms 1 or 2
options.vert_enum='C'; % use vertex enumeration variant C


[S,Sh,T,PP,PPh,DD]=bensolve(P,B,b,Y,[],c,options);

display(S);
display(Sh);
display(T);
display(PP);
display(PPh);
display(DD);

% plot upper image of (P) and lower image of (D^*)
plotresult(PP,PPh,DD,c);