function [S,Sh,T] = alg2b (P,B,b,Sh,Th,Y,~,c)
% ALG2B (subroutine of BENSOLVE) is an implementation of a modification of 
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

% store known vectors to reduce effort
Tp_known=d_hat;

% 10: T <- \emptyset (moved because of the changes)
T=[];
v=1;
while v > eps
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
      
        % 14: solve D1(w)       
        u=D1(w(t,c),P,B,b);       
        v = t(q,1) - b'*u;
        
        if v > eps
            
            % 18/19: use solution of D1(w) rather solving R^*(t^*)
            % (Webster's variant)
            s = [t(1:q-1,1);b'*u];
            
            % 20: x <- solve(P1(w(s^*)))
            x = P1(w(s,c),P,B,b);
            
            % 21: S <- S U {x}
            S = [S x]; %#ok
            break;
        else
            Tp_known=[Tp_known t]; %#ok
            
            % 15: T <- T U {(u,w(t^*))} 
            T = [T [u;w(t,c)]]; %#ok
        end
    end
end
end

function o = w(s,c)
    o=[s(1:end-1,1);1-c(1:end-1,1)'*s(1:end-1,1)];
end


function x = P1(w,P,B,b)
% solves problem (P_1(w)), cmp. [1], Section 4.1.1
%
% min  w^T Px 
% s.t. Bx>=b    
    [x,~,flag] = lpsolve(w'*P,-B,-b);
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

