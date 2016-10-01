function [S,Sh,T,PP,PPh,DD] = bensolve( P , B , b , Y , Z , c , opt )
% BENSOLVE solves linear vector optimization problems
%
% C-minimize Px
%   s.t. Bx>=b
%
% USAGE:
%
% bensolve(P,B,b)
% bensolve(P,B,b,Y ,[],c,opt)
% bensolve(P,B,b,[],Z ,c,opt)
%
% The polyhedral ordering cone C is given either by
%    Y ... matrix whose columns are the generating vectors of C
% or by
%    Z ... matrix whose columns are the generating vectors of C^+ (dual cone).
%
% The cone C is required to be polyhedral with nonempty interior and
% containing no lines. The parameter c is required to be an interior point
% of C with positive last component. This does not imply any loss of
% generality as the following equivalent problem can be considered:
% (-C)-minimize (-P)x
%      s.t. Bx>=b.
% The solution T of the dual problem (D^*) as well as the lower image of
% the (geometric) dual problem (given by DD) depend on c, compare [1].
%
% DEFAULT SETTINGS:
%
% - if no cone C is specified (by Y or Z), R^q_+ is taken by default,
%   i.e., Y=Z=eye(q)
% - if c is not specified, c=(1,...,1)^T is taken by default
% - see below for default options
%
% OUTPUT:
%
% BENSOLVE returns (compare [1] Chapter 5):
% (S,Sh) a solution to (P)
% T      a solution to (D^*)
% PP     (extreme) points of the upper image of (P)
% PPh    (extreme) directions of the upper image of (P)
% DD     (extreme) points of the lower image of (D^*)
%
% OPTIONS:
%
% - info: display information
%   0: suppress displaying information
%   1: display information (default)
%   2: display more information
%
% - dual: select primal or dual Benson algorithm
%   0: primal Benson (default)
%   1: dual Benson
%
% - alg: select variant of Benson
%   a: classical 
%   b: Webster-variant (due to an idea of Kevin Webster, Princeton)
%
% - lp_solver: select LP solver
%   0: MATLAB linprog (not recommendable for the Matlab version used)
%   1: cdd with criss-cross method, see [2] (default)
%   2: cdd with dual simplex, see [2] 
%   3: glpk (revised simplex), see [3]
%
% - eps: used in alg1 and alg2 (default 1e-07)
%
% - alg3_eps: used in alg3 (default 1e-10)
%
% - vert_enum: select the vertex enumeration method
%   A: using MATLAB convhulln, see vert1a.m and vert2a.m
%   B: using MATLAB convhulln, see vert1b.m and vert2b.m (default)
%   C: using cdd, see vert1c.m and vert2c.m
%
% - vert_enum_eps: used in vert##.m (default 1e-08)
%
% - close_vert_eps: eps to delete vertices of the upper image \P which are
%                   close to other vertices, see vert1.m and vert2.m
%                   (default 1e-12)
%
% - close_vert_eps_dual: eps to delete vertices of the upper image \D^* which
%                        are close to other vertices, see vert1.m and vert2.m
%                        (default 1e-12)
%
% SUBROUTINES
%                         bensolve.m
%                 ___________|______________________________
%                |                     |           |        |
%              alg1.m                alg2.m      alg3.m   dualcone.m
%         (alg1x.m x=a,b)      (alg2x.m x=a,b)
%                |                     |           
%                |                     |  
%             vert1.m               vert2.m 
%         (vert1x.m x=a,b,c)   (vert2x.m x=a,b,c)
%  
% further subroutines:
%  - vectorismember.m
%  - normalize.m
%  - lpsolve.m
%  - lperror.m
%  - deletemultiples.m
%  - deletezeros.m
%
% EXAMPLES:
%
% B=[2 1 1 0;1 2 0 1]';
% b=[6 6 0 0]';
% P=[1 -1; 1 1];
% [S,Sh,T,PP,PPh,DD]=bensolve(P,B,b); 
% 
% B=[2 1 1 0;1 2 0 1]';
% b=[6 6 0 0]';
% P=[1 -1; 1 1];
% c=[0 1]';
% Y=[-3 1;1 2];
% options=struct('info',2,'dual',1,'eps',1e-4);
% [S,Sh,T,PP,PPh,DD]=bensolve(P,B,b,Y,[],c,options);
%
% NOTES:
%
% - Source of the cdd and glpk mex-files:
%   http://control.ee.ethz.ch/~mpt/downloads/
% - In contrast to [1], BENSOLVE works for arbitrary pointed polyhedral
%   ordering cones C. The vector e=(1,...,1)^T used in [1] has to be
%   replaced by an interior point c of C with positive last component
%
% AUTHOR(S):
%
% written by Andreas Löhne
% URL: http://ito.mathematik.uni-halle.de/~loehne
% v 1.2; November 26, 2012
% MatLab 7.13.0.546 (R2011b)
%
% - alg1a.m is partially written by Davi Michel Valladao
% - alg1b.m is written by Anne Ehrhardt
% - two bugs in previous versions reported by Maryam Hassannasab
%
% REFERENCES:
%
% [1] Löhne, A: Vector optimization with infimum and supremum, Springer,
%     2011
% [2] Fukuda, K. CDD Homepage, http://www.ifor.math.ethz.ch/~fukuda/cdd_home/
% [3] http://www.gnu.org/s/glpk/
% [4] Barber et al.: The quickhull algorithm for convex hulls, ACM
%     Transactions of Math. Software, 1996
% [5] www.qhull.org
% [6] Heyde, F.; Löhne, A.: Geometric duality, SIAM Optimization, 2008
% 
% LICENSE:
%
% Copyright (C) 2011-2012  Andreas Löhne
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see http://www.gnu.org/licenses/.

% ---
% save PATH and set paths for LP solvers (include all subfolders)
% saved_path=path;
% addpath(genpath(pwd)); 

% count the number of LPs solved
global lp_count;
lp_count=0;

% locally used epsilons
eps1=1e-5;
eps2=1e-10;

% number of objectives q, number of variables n, number of constraints m
q = size(P,1);
[m,n] = size(B);

% set default c if necessary
if ~exist('c','var') || isempty(c)
    c=ones(q,1);
end
% last component of c must be 1
if c(q,1) > 0.01
    c=c/c(q,1);
else
    error('Last component of c must be positive.')
end

% set default options if necessary
if exist('opt','var')~=1
    opt=struct();
end
if ~isfield(opt,'info')
    opt.info=1;
end
if ~isfield(opt,'dual')
    opt.dual=0;
end
if ~isfield(opt,'variant')
    opt.variant='b';    
end
if ~isfield(opt,'lp_solver')
    opt.lp_solver=1;
end
if ~isfield(opt,'eps')
    opt.eps=1e-7;
end
if ~isfield(opt,'alg3_eps')
    opt.alg3_eps=1e-10;
end
if ~isfield(opt,'vert_enum')
    opt.vert_enum='B';
end
if ~isfield(opt,'vert_enum_eps')
    opt.vert_enum_eps=1e-8;
end
if ~isfield(opt,'close_vert_eps')
    opt.close_vert_eps=1e-12;
end
if ~isfield(opt,'close_vert_dual_eps')
    opt.close_vert_dual_eps=1e-12;
end

if opt.variant=='B'
    opt.variant='b';
end
if opt.variant=='A'
    opt.variant='a';
end
if opt.vert_enum=='a'
    opt.vert_enum='A';
end
if opt.vert_enum=='b'
    opt.vert_enum='B';
end
if opt.vert_enum=='c'
    opt.vert_enum='C';
end

global options
options=opt;

% display info
if options.info>=1
    if options.dual==0
        display('*****************************');
        display('*  Primal Benson algorithm  *');
        display('*****************************');
        display(sprintf('-- variables:     %d',n));
        display(sprintf('-- constraints:   %d',m));
        display(sprintf('-- objectives:    %d',q));
        display(sprintf('-- variant:       %c',options.variant));
        display(sprintf('-- parameter c:   %s',mat2str(c)));
        lps={'0 (MATLAB linprog)','1 (cdd criss-cross)','2 (cdd dual simplex)','3 (glpk revised simplex)'};
        display(sprintf('-- lp_solver:     %s',lps{options.lp_solver+1}));
        display(sprintf('-- eps:           %0.e',options.eps));
        display(sprintf('-- alg3_eps:      %0.e',options.alg3_eps));
        display(sprintf('-- vert_enum:     variant %c',options.vert_enum));
        display(sprintf('-- vert_enum_eps: %0.e',options.vert_enum_eps));       
    elseif options.dual==1
        display('***************************');
        display('*  Dual Benson algorithm  *');
        display('***************************');
        display(sprintf('-- variables:     %d',n));
        display(sprintf('-- constraints:   %d',m));
        display(sprintf('-- objectives:    %d',q)); 
        display(sprintf('-- variant:       %c',options.variant));
        display(sprintf('-- parameter c:   %s',mat2str(c))); 
        lps={'0 (MATLAB linprog)','1 (cdd criss-cross)','2 (cdd dual simplex)','3 (glpk revised simplex)'};
        display(sprintf('-- lp_solver:     %s',lps{options.lp_solver+1}));
        display(sprintf('-- eps:           %0.e',options.eps));
        display(sprintf('-- alg3_eps:      %0.e',options.alg3_eps));
        display(sprintf('-- vert_enum:     variant %c',options.vert_enum));        
        display(sprintf('-- vert_enum_eps: %0.e',options.vert_enum_eps));        
    end
end
if options.info>=2
    display('calling BENSOLVE ...');
end

% Compute generating vectors of C or generating vectors of the dual
% cone C^+ depending on what is given: 

% C = R^q_+ is the default ordering cone
if (~exist('Y','var') && ~exist('Z','var')) || (isempty(Y) && isempty(Z))
    Y=eye(q);
    Z=eye(q);
% if only Z is given, compute Y    
elseif isempty(Y)  
    Z=deletezeros(Z,1e-10);
    Z=Z./(ones(q,1)*c'*Z); % Note: c^T*z=1 for the columns of Z
    Y=dualcone(Z);
% if only Y is given (or both Y and Z), compute Z     
else 
    Y=deletezeros(Y,1e-10);
    Y=normalize(Y);
    Z=dualcone(Y);
    Z=Z./(ones(q,1)*c'*Z); % Note: c^T*z=1 for the columns of Z
end  

% Compute eta using Algorithm 3, see [1] Chapter 5
if options.info>=1
    display('  alg3: computing interior point of D^* ...');
end    
p = size(Y,2);
r=alg3([zeros(q-1,m),eye(q-1),zeros(q-1,1)],...
    [-eye(m),zeros(m,q);zeros(p,m),-Y'],zeros(m+p,1),...
    [[B',-P'];[zeros(1,m),c']],[zeros(n,1);1]);
eta=[r;1-c(1:q-1,1)'*r];

% solve the homogeneous problem, compare Section 5.4 in [1]
if options.info>=1
    display('  done.');
    display('  alg1: solving homogeneous problem ...');
end
if options.dual==0
    if options.variant=='a'
       [S_eta,~,T_eta]=alg1a(P,[B;-eta'*P],[zeros(m,1);-1],[],...
           [zeros(m+1,size(Z,2));Z],Y,Z,c); % classical variant
    elseif options.variant=='b' 
       [S_eta,~,T_eta]=alg1b(P,[B;-eta'*P],[zeros(m,1);-1],[],...
           [zeros(m+1,size(Z,2));Z],Y,Z,c); % Webster variant
    else error('Variant choesn for Algorithm 1 does not exist.');
    end   
elseif options.dual==1
    if options.variant=='a'
       [S_eta,~,T_eta]=alg2a(P,[B;-eta'*P],[zeros(m,1);-1],[],...
           [zeros(m+1,size(Z,2));Z],Y,Z,c); % classical variant
    elseif options.variant=='b'
       [S_eta,~,T_eta]=alg2b(P,[B;-eta'*P],[zeros(m,1);-1],[],...
           [zeros(m+1,size(Z,2));Z],Y,Z,c); % Webster variant
    else
        error('Variant choesn for Algorithm 2 does not exist.');
    end    
end

% construct Sh, see Th. 5.27 in [1]
Sh=[];
for i=1:size(S_eta,2)
    x=S_eta(:,i);    
    if norm(P*x) > eps1        
        Sh=[Sh,x/norm(x)];%#ok
    end
end

% construct Th, see Th. 5.28 in [1]
Th=[];
for i=1:size(T_eta,2)
    v=T_eta(:,i);
    if abs(v(m+1,1)) < eps1
        Th =[Th,[v(1:m,1);v(m+2:m+q+1,1)]]; %#ok
    end
end

% solve the inhomogeneous problem
if options.info>=1
    display('  done.');
    display('  alg1: solving inhomogeneous problem ...');
end

if options.dual==0 % use primal Benson
    if options.variant=='a'
       [S,Sh,T]=alg1a(P,B,b,Sh,Th,Y,Z,c); % classical variant
    elseif options.variant=='b'
       [S,Sh,T]=alg1b(P,B,b,Sh,Th,Y,Z,c); % Webster variant
    else
        error('Variant chosen for Algorithm 1 does not exist.');
    end        
elseif options.dual==1 % use dual Benson
    if options.variant=='a'
       [S,Sh,T]=alg2a(P,B,b,Sh,Th,Y,Z,c); % classical variant
    elseif options.variant=='b'
       [S,Sh,T]=alg2b(P,B,b,Sh,Th,Y,Z,c); % Webster variant
    else
        error('Variant chosen for Algorithm 2 does not exist.');
    end 
end

if options.info>=1
    display('  done.');    
end

% return the upper image (PP,PPh) of (P) and the lower image DD of (D^*)
PP=deletemultiples(P*S,eps2);
if ~isempty(Sh)
    PPh=deletemultiples([P*Sh,Y],eps2);
else
    PPh=deletemultiples(Y,eps2);
end
DD=deletemultiples([T(m+1:m+q-1,:); b' * T(1:m,:)],eps2);

if options.info>=2 
    display(sprintf('Vector optimization problem solved. (%d LPs).',lp_count));
end    
lp_count=0;

% restore PATH
% path(saved_path);

end



