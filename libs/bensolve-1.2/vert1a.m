function [ Tp ] = vert1a( Td, c )
% Subroutine for Algorithm 1 in chapter 5 of [1] (primal Benson algorithm)
% transforms T_d in T_p (line 06 of Algorithm 1, vertex enumeration)
%
% VARIANT A: The matlab function convhulln is used, which is based on QHULL,
% compare [4,5]. The transformations are based on "geometric duality" [6,1]
% 
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m

% tolerance
global options
eps =options.vert_enum_eps; 

% tolerance used in vectorismember()
eps2=1e-10;  

% Td consits of p q-dimensional column vectors
[q,p]=size(Td);

% The convex hull conv(Td) of the column vectors of Td may have an empty
% interior, whereas the interior of conv(Td)-K for K={t*e^q| t<=0} is
% nonempty. Since the matlab function convhulln only works with bounded
% polyhedra, we use conv(Td) instead of conv(Td)-K. Additional facets will
% be ruled out below. To avoid an empty interior of conv(Td), we append an
% appropriate vector to Td
Td=[Td,[sum(Td(1:q-1,:),2)/p;min(Td(q,:)-3.14)]];
p=p+1;

% Convhulln returns a list of indices that describes which column vectors 
% of Td define a facet of conv(Td)
Ix=convhulln(Td');

% r is the number of "facets" of conv(Td) each of which is defined by
% exactly q vertices, note that q is not changed here
[r,q]=size(Ix);
Tp=[];

% Each facet of conv(Td) is tested for being K-maximal, if so it
% corresponds by geometric duality [4] to a vertex of \P, which is stored in
% Tp. 
for j=1:r 
    
    % Those q vectors of Td (of dimension q) that define the j-th "facet" of
    % conv(Td) are stored in Z.
    Z=Td(:,Ix(j,:))';  

    % Compute a normal vector of the j-th facet.
    nv = null(Z(1:q-1,:)-ones(q-1,1)*Z(q,:));

    % Only non-vertical facets are used; the second condition tests if we
    % have a facet (i.e. the dimension is q-1 and thus the dimension of
    % the nullspace is 1, faces with dimension less than q-1 occured in an
    % example with q=4); the first condition ensures that the facet is
    % not vertical          
    if(abs(nv(q))>eps && size(nv,2)==1)        
        
        % Scale the normal vector nv of the j-th facet such that the last
        % component is -1.
        nv=-nv./nv(q);
        
        % Compute the righthand side d of the supporting hyperplane
        % {y^*|nv^T y^*=d} of conv(Td) that contains the j-th facet.
        d=Z(1,:)*nv;
        
        % Test whether the j-th facet is K-maximal or not: try to find a vector
        % y^* of Td such that nv^T y^* > d, once such a vector found, we know
        % the j-th facet is K-maximal. Then, the equality of the 
        % corresponding supporting hyperplane is used to compute the
        % elements y of Tp. The relationsship is given as nv=w^*(y) and
        % d=y_q, where w^*(y) := (y_1- c_1*y_q, ... y_{q-1}-c_{q-1}*y_q,-1),
        % compare (4.16) in [1] for the case that c=(1,...,1)^T
        for i=1:p
            if (nv' * Td(:,i) - d > eps2)
                e=[nv(1:q-1,1)-d*c(1:q-1,1);-d];  
                if(~vectorismember(e,Tp,eps2))
                    Tp=[Tp,e]; %#ok

                end                               
                break
            end
        end        
    end
end
end
