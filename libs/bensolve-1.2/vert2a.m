function [Tp] = vert2a(Td,Td_hat,Y,c)
% Subroutine of Algorithm 2 in Chapter 5 of [1] (dual Benson algorithm)
% transforms (T_d,\hat T_d) in T_p (line 09 of Algorithm 2, vertex enumeration)
%
% VARIANT A: The matlab function convhulln is used, which is based on QHULL,
% compare [4,5]. The transformations are based on "geometric duality" [6,1]
% 
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m

% tolerance
global options
eps=options.vert_enum_eps;

% tolerance used in vectorismember()
eps2=1e-10;

% Append the generating vectors of C to Td_hat (they are by definition not
% contained)
Td_hat=[Td_hat Y];

% Td consits of p q-dimensional column vectors
[q,p]=size(Td);
% Td_hat consits of ph q-dimensional column vectors
ph=size(Td_hat,2);

% Compute an interior point d_hat of conv(Tp)-K, where K=R_+.{e^q}, in
% two steps: 
% Step 1: Compute the first q-1 components of d_hat. Consider the
% (q-1)-dimensional polytop given by the inequalities 
% \varphi^h(y,y^*)>=0, where y are the column vectors of Td_hat
% (see Chapter 4.6 in [1]):
d_hat=IntPoint(-(Td_hat(1:q-1,:)-c(1:q-1,1)*Td_hat(q,:))',Td_hat(q,:)');

% Step 2: Compute the q-th component of d_hat directly by ensuring that 
% w^*(y)^T d_hat > -y_q for all vectors y \in Td 
d_hat(q,1)=min((d_hat'*(Td(1:q-1,:)-c(1:q-1,1)*Td(q,:))+Td(q,:)),[],2)-1;

% The idea of the next three steps is as follows: From Td and Td_hat we
% obtain by geometric duality an inequality representation of the dual
% polyhedron D^*. But CONVHULLN needs vectors as input. D^* is regarded to
% be a polytope with a vertex at infinity and the dual polytope is
% computed. To this end, the inequalities are transformed such that its
% r.h.s is 1 and the ordering is <= (steps 1 and 2). The coefficient
% vectors are the vertices of the dual polytope.

% Step 1: V consists of the vectors w^*(y)^T /(-y_q - w^*(y)^T d_hat))
% where y \in Td (see [1] Chapter 4.5 for notation) 
V=[Td(1:q-1,:)-c(1:q-1,1)*Td(q,:);-ones(1,p)]*...
  diag(-1./([Td(1:q-1,:)-c(1:q-1,1)*Td(q,:);-ones(1,p)]'*d_hat+Td(q,:)'));

% Step 2: Vh consists of the vectors w^{*h}(y)^T /(-y_q - w^{*h}(y)^T d_hat))
% where y \in Td_hat (w^{*h}(y)^T is definined like w^*(y)^T but the
% last component is zero)
Vh=[Td_hat(1:q-1,:)-c(1:q-1,1)*Td_hat(q,:);zeros(1,ph)]*...
   diag(-1./([Td_hat(1:q-1,:)-c(1:q-1,1)*Td_hat(q,:);zeros(1,ph)]'...
   *d_hat+Td_hat(q,:)'));

% Step 3: U contains the vertices of the dual polytope of D^*, which is
% obtained from D^* by the usual duality (polarity) map for polytopes.
% Although D^* is not bounded, this transformation can be used. The only
% differrence is that the dual polytope contains zero in one of its facets.
% This facets corresponds to a vertex of D^* at infinity and has to be
% ruled out later.
U=[V,Vh]';

% compute indices of points in U that define "facets" of the dual
% polytope of D^*, such facets correspond to the vertices of D^*
Ix=convhulln(U);
                
[r,s]=size(Ix); Tp=[];
for j=1:r 
    % points in U that define the j-th facet of D^** are stored in B
    B=U(Ix(j,:),:); 
    
    % in particular, the facet of the dual polytope of D^* that contains
    % the origin is ruled out here
    if(abs(det(B))>eps)
        % compute the coefficients to the inequalities describing the
        % facets of the dual polytope of D^*
        c=B\ones(s,1); 
        
        % the coefficients correspond (up to a shift by d_hat) to the
        % vertices of D^*
        e=c+d_hat;
        if(~vectorismember(e,Tp,eps2))
            Tp=[Tp,e]; %#ok
        end 
    end                
end
end

% compute a relative interior point of a polytope given by inequalities
function s = IntPoint(B , b )                           
            [x,~,flag] =lpsolve([zeros(size(B,2),1);1],...
                [B';-ones(1,size(B,1))]',b);
            lperror(flag);    
            s=x(1:end-1);            
end
        






