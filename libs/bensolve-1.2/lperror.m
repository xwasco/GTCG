function lperror(flag)
% LPERROR produces an error message for different LP solvers
%
% LICENSE: This file is part of BENSOLVE, see bensolve.m

global options

%for MATLAB linprog
if options.lp_solver==0
    if flag == 0
        error('linprog error: Maximum number of iterations reached.')
    elseif flag == -2
        error('linprog error: No feasible point found.')
    elseif flag == -3
        error('linprog error: Problem is unbounded.')
    elseif flag == -4
        error('linprog error: NaN value encountered during execution of algorithm.')
    elseif flag == -5
        error('linprog error: Both primal and dual problems are infeasible.')
    elseif flag == -7
        error('linprog error: Magnitude of search direction became too small; no further progress can be made. The problem is ill-posed or badly conditioned.')
    end

% for cdd LP solver
elseif options.lp_solver==1 || options.lp_solver==2
    if flag == 0
        error('cdd error: dd_LPSundecided')
    elseif flag == 2
        error('cdd error: dd_Incosistent')
    elseif flag == 3
        error('cdd error: dd_DualIncosistent')
    elseif flag == 4
        error('cdd error: dd_StrucIncosistent')
    elseif flag == 5
        error('cdd error: dd_StrucDualIncosistent')
    elseif flag == 6
        error('cdd error: dd_Unbounded')
    elseif flag == 7
        error('cdd error: dd_DualUnbounded')
    end
% for glpk LP solver  
elseif options.lp_solver==3
    if flag == 1
        error('glpk error: solution is undefined.')
    elseif flag == 2
        error('glpk error: solution is feasible.')
    elseif flag == 3
        error('glpk error: solution is infeasible.')
    elseif flag == 4
        error('glpk error: no feasible solution exists.')
    elseif flag == 6
        error('glpk error: solution is unbounded.')
    elseif flag == 101
        error('glpk error: invalid basis.')
    elseif flag == 102
        error('glpk error: singular matrix.')
    elseif flag == 103
        error('glpk error: ill-conditioned matrix.')
    elseif flag == 104
        error('glpk error: invalid bounds.')    
    elseif flag == 105
        error('glpk error: solver failed.')
    elseif flag == 106
        error('glpk error: objective lower limit reached.')
    elseif flag == 107
        error('glpk error: objective upper limit reached.')    
    elseif flag == 108
        error('glpk error: iteration limit exceeded.')
    elseif flag == 109
        error('glpk error: time limit exceeded.')
    elseif flag == 110
        error('glpk error: no primal feasible solution.')
    end   
end
end

