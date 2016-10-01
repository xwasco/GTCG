function [ Tp ] = vert2c( Td , Td_hat, d_hat ,Y , c )
% Subroutine of Algorithm 2 in Chapter 5 of [1] (primal Benson Algorithm)
% transforms T_d in T_p (line 09 of Algorithm 2, vertex enumeration)
% 
% VARIANT C : The dual representation of the polyhedron is transformed
% into a polytope. cdd [2] is used to compute a dual polytope. The dual
% polytope is transferred into primal representation of the polyhedron.
%
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m

%%% Transform polyhedron into a polytope with the same facial structure

% inequalities of non-verticel hyperplanes
[q,p1]=size(Td);
C1=[Td(1:q-1,:)-c(1:q-1,1)*Td(q,:);-ones(1,p1)]';
d1=-Td(q,:)';

% Append the generating vectors of C to Td_hat (they are by definition not
% contained)
Td_hat=[Td_hat Y];

% inequalities of vertical hyperplanes
[q,p2]=size(Td_hat);
C2=[Td_hat(1:q-1,:)-c(1:q-1,1)*Td_hat(q,:);zeros(1,p2)]';
d2=-Td_hat(q,:)';

C=[C1;C2];
d=[d1;d2];
p=p1+p2;

M=C'*diag(1./(d-C*d_hat));

xi=sum(M,2)/(1+size(M,2));

db=M-xi*ones(1,p);

%%% Compute dual polytope and rule out the facet that corresponds to an
%%% extreme direction 

global options
eps=options.vert_enum_eps;

H=struct('A',db','B',ones(p,1));
V=cddmex('extreme',H);
pb1=V.V';

pb=[];
for i=1:size(pb1,2)
    if 1 + pb1(:,i)'*xi > eps
        pb=[pb,pb1(:,i)]; %#ok
    end
end    

%%% Reverse transformation

p=size(pb,2);
Tp=d_hat*ones(1,p)+pb*diag(1./(ones(p,1) + pb'*xi));

end


