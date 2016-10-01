% simplest example

B=[2 1 1 0;1 2 0 1]';
b=[6 6 0 0]';
P=[1 -1; 1 1];

[S,Sh,T,PP,PPh,DD]=bensolve(P,B,b);

display(S);
display(Sh);
display(T);
display(PP);
display(PPh);
display(DD);
