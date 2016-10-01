function [ Tp ] = vert1b( Td , p_hat, c )
% Subroutine for Algorithm 1 in chapter 5 of [1] (primal Benson algorithm)
% transforms T_d in T_p (line 06 of Algorithm 1, vertex enumeration)
%
% VARIANT B: The matlab function convhulln is used, which is based on QHULL,
% compare [4,5]. The Transformations are the same as in Variant C
%
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m


%%% Transform polyhedron to polytope

[q,p]=size(Td);

C=[Td(1:q-1,:);ones(1,p)-c(1:q-1,1)'*Td(1:q-1,:)]';
d=Td(q,:)';

M=C'*diag(1./(d-C*p_hat));

xi=sum(M,2)/(1+size(M,2));

db=[M-xi*ones(1,p), -xi];

%%% Compute dual polytope and rule out faces containing the (p+1)-th vertex 

% tolerance
global options
eps=options.vert_enum_eps;

% tolerance used in vectorismember()
eps2=1e-10;  

Ix=convhulln(db');

pb=[];

for j=1:size(Ix,1) 
    
    % this condition could be omitted as there is a second condition ruling
    % out faces containing the (p+1)-th vertex
    if ~ismember(p+1,Ix(j,:))
        
        % Those q vectors of pd (of dimension q) that define the j-th "facet"
        % of conv(pd) are stored in Z.        
        Z=db(:,Ix(j,:))';     

        % Compute a normal vector of the j-th facet.
        c = null(Z(1:q-1,:)-ones(q-1,1)*Z(q,:)); 

        % this condition ensures that we have a facet (and not a lower
        % dimensional face) 
        if size(c,2)==1 
            % Compute the righthand side d of the supporting hyperplane
            % {y^*|c^T y^*=d} of conv(vp) that contains the j-th facet.
            d=Z(1,:)*c; 
            e=c/d;  
            % the second condition rules out faces containing vertex p+1
            % in case it is not indicated by indices from convhulln
            if(~vectorismember(e,pb,eps2)&&(1 + e'*xi)>eps)
                pb=[pb,e]; %#ok
            end            
        end
    end

end

%%% Reverse transformation

p=size(pb,2);
Tp=p_hat*ones(1,p)+pb*diag(1./(ones(p,1) + pb'*xi));

end

