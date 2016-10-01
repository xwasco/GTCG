function [S,Sh,T] = alg1b(P,B,b,Sh,Th,Y,Z,c)
% ALG1B (subroutine of BENSOLVE) is an implementation of a modification 
% (Webster variant) of Algorithm 1 in [1]
%
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m
 
global options
eps=options.eps;

q = size(P,1);
m = size(B,1);

% 02: p_hat <- P(solve(P1(0)) + c
p_hat = P*P1(zeros(q,1),P,B,b) + c;

% 03: T <- {(solve(D1(w)),w)|(u,w) \in Th}
for i = 1:size(Th,2)
    w = Th(m+1:m+q,i);   
    T(1:m+q,i) = [D1(w,P,B,b);w];  %#ok
end

z=1;

% store known vectors to reduce the effort
Tp_known=[];

% 07: S <- emptyset 
S = [];

while z>eps
    % 05: Td <- {D*(u,w)|(u,w) \in T}
    Td=[T(m+1:m+q-1,:); b' * T(1:m,:)];
    
    % 06: Tp -> vert(Td)
    Tp = vert1(Td,p_hat,c);
    
    % select new vectors in Tp
    % Tp_unknown=setdiff(Tp',Tp_known','rows')';    
    
    Tp_unknown=[];
    for i=1:size(Tp,2)
         if ~vectorismember(Tp(:,i),Tp_known,1e-13)
             Tp_unknown=[Tp_unknown,Tp(:,i)]; %#ok
         end
    end  
    
    if isempty(Tp_unknown)    % new in version 1.2   
        break;                % thanks to Maryam Hassannasab for report
    end                       % 
    
    for i = 1:size(Tp_unknown,2)
        % 10: t<-Tp[i]
        t = Tp_unknown(1:q,size(Tp_unknown,2)-i+1);         
       
        [x,z]=P2(t,P,B,b,Z,c);      
        
        % 13: replaced by the test z > 0
        if z>eps
            % 16: s <- t + z*c
            s = t + z*c;
            
            % 17: (u,w) <- solve(D2(s))
            [u,w] = D2(s,P,B,b,Y,c);        
            
            % 18: T <- T U {(u,w)}
            T = [T [u;w]];  %#ok
            break;
        else
            Tp_known=[Tp_known t]; %#ok
            
            % 12: S <- S U {x} 
            S = [S x];  %#ok
        end
    end
end
end

function x = P1(w,P,B,b)
% solves problem (P_1(w)), cmp. [1], Section 4.1.1
%
% min  w^T Px 
% s.t. Bx>=b    
    [x,~,flag] = lpsolve((w'*P)',-B,-b);
    lperror(flag);   
end

function u = D1(w,P,B,b)
% solves problem (D_1(w)), cmp. [1], Section 4.1.1
%
% max  b^T u
% s.t. B^T u  = P^T w
%          u >= 0   
    [u,~,flag] = lpsolve(-b,[],[],B',P'*w,zeros(length(b),1));
    lperror(flag);
end

function [x,z] = P2(t,P,B,b,Z,c)
% min  z
% s.t. B*x >= b
%      Z^T*P*x <= Z^T*(y + z*c)
    
    n = size(P,2);
    m = size(B,1);    
    [aux,~,flag] = lpsolve([zeros(n,1); 1],[-B zeros(m,1); Z'*P -Z'*c],[-b;Z'*t]);
    x = aux(1:n,1);
    z = aux(n+1,1);
    lperror(flag); 
end

function [u,w] = D2(y,P,B,b,Y,c)
% solves problem (D_2(y)), cmp. [1], Section 4.1.1
% Y is added to handle arbitrary ordering cones
%
% max  b^T u - y^T w
% s.t. B^T u - P^T w  = 0
%              c^T w  = 1
%                  u >= 0
%               Y' w >= 0 
    
    [q,n] = size(P);
    m = size(B,1);
    p = size(Y,2);    
    [aux,~,flag] = lpsolve(-[b ; -y],[-eye(m),zeros(m,q);zeros(p,m),-Y'],zeros(m+p,1),[B' -P'; zeros(1,m) c'],[zeros(n,1);1]);
    u = aux(1:m,1);
    w = aux(m+1:m+q,1);    
    lperror(flag);
end











