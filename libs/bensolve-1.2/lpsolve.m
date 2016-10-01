function [ x,fval,exitflag ] = lpsolve(f,A,b,Aeq,beq,lb,ub)
% LPSOLVE chooses LP solver depending on global variable options.lp_solver
% 
% lpsolver==0: MATLAB linprog
% lpsolver==1: cdd criss cross method (cdd mex-file required)
% lpsolver==2: cdd dual simplex (cdd mex-file rquired)
% lpsolver==3: glpk (glpk mex-file required)
%
% LICENSE: This file is part of BENSOLVE, see bensolve.m

if ~exist('A','var')
    A=[];
end
if ~exist('b','var')
    b=[];
end
if ~exist('Aeq','var')
    Aeq=[];
end
if ~exist('beq','var')
    beq=[];
end
if ~exist('lb','var')
    lb=[];
end
if ~exist('ub','var')
    ub=[];
end

% transpose f to row vector
if size(f,1) ~= 1 
        f=f';
end 

global options

% use MATLAB linprog 
if options.lp_solver==0    
    [x,fval,exitflag] = linprog(f,A,b,Aeq,beq,lb);
% use cdd LP solver with criss cross method
elseif options.lp_solver==1        
    IN=struct('obj',f,'A',[Aeq;A;-eye(size(lb,1));eye(size(ub,1))],'B',[beq;b;-lb;ub],'lin',(1:size(Aeq,1)));
    OUT=cddmex('solve_lp',IN);
    x=OUT.xopt;
    fval=OUT.objlp;
    exitflag=OUT.how;
% use cdd LP solver with dual simplex method
elseif options.lp_solver==2       
    IN=struct('obj',f,'A',[Aeq;A;-eye(size(lb,1));eye(size(ub,1))],'B',[beq;b;-lb;ub],'lin',(1:size(Aeq,1)));
    OUT=cddmex('solve_lp_DS',IN);
    x=OUT.xopt;
    fval=OUT.objlp;
    exitflag=OUT.how;
% use GLPK    
elseif options.lp_solver==3  
    m1=size(A,1);
    m2=size(Aeq,1);
    n=size(f,2);
    [x,fval,exitflag]=glpk(f',[A;Aeq],[b;beq],lb,ub, char([zeros(1,m1)+'U' zeros(1,m2)+'S']),char(zeros(1,n)+'C'),1);     
else
    error('No valid LP solver choosen.')
end

% count the number of LPs solved by Benson's algorithm
global lp_count
lp_count=lp_count+1;

end

