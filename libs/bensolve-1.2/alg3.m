function [ w ] = alg3(P,  B , b , Beq, beq)
% ALG3 (subroutine of BENSOLVE) is an implementation of Algorithm 3 in [1]
%
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m

global options
eps=options.alg3_eps;

w=[];
k=size(P,1);
% line 03
z=solveQ(zeros(k,1),P,B,b,Beq,beq);
% line 04
v=zeros(k,1);
c=[];
for i=1:k
    % line 7  (in contrast to [1] the smallest index is 1 rather than 0)  
    if ~isempty(c) 
        c(:,i)=v(:,i)-sum(c(:,2:end-1)*diag(1./sum(c(:,2:end-1).^2))*diag(v(:,i)'*c(:,2:end-1)),2); %#ok
    else
       c=v(:,i); 
    end    
    % line 8
    normal=null(c');
    c=[c,normal(:,1)]; %#ok   
    % line 9    
    v= [v,solveQ(c(:,i+1),P,B,b,Beq,beq)-z]; %#ok       
    % lines 10 and 11
    if abs(v(:,i+1)'*c(:,i+1))<eps 
        v(:,i+1)= solveQ(-c(:,i+1),P,B,b,Beq,beq)-z;        
    end    
    % lines 12 to 16 
    if abs(v(:,i+1)'*c(:,i+1))<eps         
        return;
    end
end
%line 18
w=z+ 1/(k+1).*sum(v,2); 
end

function [y] = solveQ(c,P,B,b,Beq,beq)
% min c^ x
% s.t. B  *x <= b 
%      Beq*x == beq
    [x,~,flag] = lpsolve(P'*c,B,b,Beq,beq);
    lperror(flag);    
    y=P*x;
end


