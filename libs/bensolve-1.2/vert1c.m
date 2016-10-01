function [ Tp ] = vert1c( Td , p_hat , c )
% Subroutine of Algorithm 1 in chapter 5 of [1] (primal Benson algorithm)
% transforms T_d in T_p (line 06 of Algorithm 1, vertex enumeration)
% 
% VARIANT C : The dual representation of the polyhedron is transformed
% into a polytope. cdd [2] is used to compute a dual polytope. The dual
% polytope is transformed into primal representation of the polyhedron.
%
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m

%%% Transform polyhedron into a polytope with the same facial structure

[q,p]=size(Td);

C=[Td(1:q-1,:);ones(1,p)-c(1:q-1,1)'*Td(1:q-1,:)]';
d=Td(q,:)';

M=C'*diag(1./(d-C*p_hat));

xi=sum(M,2)/(1+size(M,2));

db=[M-xi*ones(1,p), -xi];

%%% Compute dual polytope and rule out faces containing the (p+1)-th vertex 

global options
eps =options.vert_enum_eps;

H=struct('A',db','B',ones(p+1,1));
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
Tp=p_hat*ones(1,p)+pb*diag(1./(ones(p,1) + pb'*xi));

end


