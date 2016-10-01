% Example with q=3 and 4 generating vectors of C

B=[eye(3);ones(1,3);1 2 2;2 2 1;2 1 2];
b=[0 0 0 1 3/2 3/2 3/2]';
P=[1 1 0;0 1 1;1 0 1]';

% generating vectors of ordering cone C
Y=[1 0 0 ; 0 1 0 ; -1 0 2 ; 0 -1 2]';

clear options;
options.vert_enum='C';
[S,Sh,T,PP,PPh,DD]=bensolve(P,B,b,Y,[],[],options);

display(PP);
display(PPh);
display(DD);

plotresult(PP,PPh,DD);


