function [ Z ] = dualcone ( Y )
% DUALCONE computes the generating vectors Z of the dual cone of a
% polyhedral cone C which is given by its generating vectors Y
%
% LICENSE: This file is part of BENSOLVE, see bensolve.m

global options

% choose an LP solver
if ~isfield(options,'lp_solver')
    options.lp_solver=1;
end

% locally used epsilon
eps=1e-10; 

Y=deletezeros(Y,1e-10);

[q,p]=size(Y);

if q==1 % the one-dimensional case
    c1=min(Y);
    c2=max(Y);
    if c2 * c1 < 0
        Z=0;
    else 
        Z=c1/norm(c1);
    end
else % the multi-dimensional case

    % compute interior point of the dual cone
    eta=IntPoint_ub(-Y',zeros(p,1));

    % compute vertices of a basis of C
    M=Y*diag(1./(Y'*eta));

    % compute an interior point of conv (M \cup {0})
    xi=sum(M,2)/(1+size(M,2));

    % change coordinates so that that 0 belonts to interior of conv (M \cup {0})
    db=[M-xi*ones(1,p), -xi];

    % Compute dual polytope and take only faces containing the (p+1)-th vertex 
    Ix=convhulln(db');

    Z=[];

    for j=1:size(Ix,1) 

        if ismember(p+1,Ix(j,:))

            % Those q vectors of pd (of dimension q) that define the j-th "facet"
            % of conv(pd) are stored in Z.        
            X=db(:,Ix(j,:))';     

            % Compute a normal vector of the j-th facet.
            c = null(X(1:q-1,:)-ones(q-1,1)*X(q,:));
            
            if X(q,:)*c > 0
                c=-c;
            end

            % this condition ensures that we have a facet (and not a lower
            % dimensional face) 
            if size(c,2)==1            
                if ~vectorismember(c,Z,eps)
                    Z=[Z,c/norm(c)]; %#ok
                end            
            end
        end
    end
end
end

% computes a relative interior point of a unbounded polyhedron given 
% by inequalities Bx<=b
function s = IntPoint_ub(B , b )                       
            [x,~,flag] =lpsolve([zeros(size(B,2),1);1],...
                [[B';-ones(1,size(B,1))]';[zeros(1,size(B,2)) -1]],[b;1]);
            lperror(flag);    
            s=x(1:end-1);            
end



