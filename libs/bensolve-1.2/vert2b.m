function [ Tp ] = vert2b( Td , Td_hat, d_hat ,Y , c )
% Subroutine of Algorithm 2 in Chapter 5 of [1] (primal Benson algorithm)
% transforms T_d in T_p (line 09 of Algorithm 2, vertex enumeration)
% 
% VARIANT B: The matlab function convhulln is used, which is based on QHULL,
% compare [2,3]. The transformations are the same as in Variant C
%
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m

%%% Transform polyhedron into a polytope with the same facial structure

% inequalities of non-vertical hyperplanes
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

% tolerance
global options
eps=options.vert_enum_eps;

% tolerance used in vectorismember()
eps2=1e-10;  

Ix=convhulln(db');

pb=[];

for j=1:size(Ix,1)         
    % Those q vectors of db (of dimension q) that define the j-th "facet"
    % of conv(db) are stored in Z.        
    Z=db(:,Ix(j,:))';     

    % Compute a normal vector of the j-th facet.
    c = null(Z(1:q-1,:)-ones(q-1,1)*Z(q,:)); 

    % this condition ensures that we have a facet (and not a lower
    % dimensional face) 
    if size(c,2)==1 
        % Compute the righthand side d of the supporting hyperplane
        % {y^*|c^T y^*=d} of conv(db) that contains the j-th facet.
        d=Z(1,:)*c; 
        e=c/d;  
        % this condition rules out faces corresponding to extreme directions
        if(~vectorismember(e,pb,eps2)&&(1 + e'*xi)>eps)
            pb=[pb,e]; %#ok
        end            
    end
end   

%%% Reverse transformation

p=size(pb,2);
Tp=d_hat*ones(1,p)+pb*diag(1./(ones(p,1) + pb'*xi));

end


