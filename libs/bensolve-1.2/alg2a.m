function [S,Sh,T] = alg2a(P,B,b,Sh,Th,Y,~,c)
% ALG2A (subroutine of BENSOLVE) is an implementation of a modification of 
% Algorithm 2 in [1]
%
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m

global options
eps=options.eps;

q = size(P,1);
m = size(B,1);
p=size(Th,2);

% 02: eta: mean of the w-components of Th
d_hat=1/p*sum(Th(m+1:m+q,:),2);

% 03: x <- solve(P_1(eta))
x=P1(d_hat,P,B,b);

% 04: d_hat <- {eta_1,...,eta_{q-1},eta^Px-1}
d_hat(q,1) = d_hat'*P*x - .1;

% 05: S <- {x}
S = x;

% 06: Td_hat <- {Px | x \in Sh}
if ~isempty(Sh)
    Td_hat=P*Sh;
else
    Td_hat=[];
end

alpha=0;
% store known vectors to reduce effort
Tp_known=d_hat;

% 10: T <- \emptyset
T=[];

while alpha < 1-eps
    % 08: Td <- {Px|x \in S}
    Td=P*S;    
      
    % 09: Tp <- vert(Td,Td_hat)
    Tp = vert2(Td,Td_hat,d_hat,Y,c);
    
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
        % 13: t^*<-Tp[i]
        t = Tp_unknown(1:q,size(Tp_unknown,2)-i+1);           
        
        % 14: D_1(w) replaced by R^*(t^*)
        [u,~,alpha]=R_s(t,d_hat,P,B,b,Y,c);
        
        % 16: replaced by the test alpha < 1
        if alpha < 1-eps
            % 18: (u,w,alpha) <- solve(R^*(t^*))
            % already solved above
            
            % 19: s^* <- alpha*t + (1-alpha)*p_hat
            s = alpha*t + (1-alpha)*d_hat;
            
            % 20: x <- solve(P1(w(s^*)))
            x = P1(ww(s,c),P,B,b);
            
            % 21: S <- S U {x}
            S = [S x]; %#ok
            break;
        else
            Tp_known=[Tp_known t]; %#ok
            
            % 15: T <- T U {(u,w(t^*))} (moved because of changes)
                       
            T = [T [u;ww(t,c)]]; %#ok
        end
    end
end
end

function o = ww(s,c)
o=[s(1:end-1,1);1-c(1:end-1,1)'*s(1:end-1,1)];
end


function x = P1(w,P,B,b)
% solves problem (P_1(w)), cmp. [1], Section 4.1.1
%
% min  w^T Px 
% s.t. Bx>=b    
    [x,~,flag] = lpsolve((w'*P)',-B,-b);
    lperror(flag);   
end

function [u,w,alpha] = R_s(t,d_hat,P,B,b,Y,c)
% solves problem (R^*(t^*)), cmp. [1], Section 5.2
% Y is added to handle arbitrary ordering cones
%
% max  alpha_s
% s.t.                B^T u - P^T w  = 0
%                             c^T w  = 1
%                                 u >= 0
%                              Y' w >= 0
%      alpha^* t^* + (1-alpha) d_hat = D^*(u,w)
    
    [q,n] = size(P);
    m = size(B,1);
    p = size(Y,2);   
       
    [aux,~,flag] = lpsolve([zeros(m+q,1);-1],...
        [-eye(m),zeros(m,q+1);zeros(p,m),-Y',zeros(p,1);zeros(1,m+q+1)],...
        zeros(m+p+1,1),[B' -P' zeros(n,1);...
        zeros(1,m) c' 0;...
        [zeros(q-1,m) eye(q-1) zeros(q-1,1); b' zeros(1,q)] d_hat-t],...
        [zeros(n,1);1;d_hat]);
    u = aux(1:m,1);
    w = aux(m+1:m+q,1);
    alpha= aux(m+q+1,1);
    lperror(flag); 
end
