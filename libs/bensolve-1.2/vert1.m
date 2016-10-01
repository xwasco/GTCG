function [ Tp ] = vert1(Td,p_hat,c)
% VERT1 chooses a vertex enumeration method for the primal Benson algorithm
%
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m

global options

% delete vertices of the outer approximation Td of the lower image of (D^*)
% being very close to each other (for stability reasons) 
Td=deletemultiples(Td,options.close_vert_dual_eps);

if size(Td,1)==1 % the scalar case
    Tp=max(Td);
else    
    % Choose the variant of the vert1 subroutine 

    % based on convhulln and dual polyhedra
    if options.vert_enum=='A'
        Tp=vert1a(Td,c);

    % based on convhulln and polyhedron-polytope transformation
    elseif options.vert_enum=='B'
        Tp=vert1b(Td,p_hat,c);        

    % based on cdd and polyhedron-polytope transformation
    elseif options.vert_enum=='C'
        Tp=vert1c(Td,p_hat,c);    
    else
        error('Chosen variant for vertex enumeration does not exist.');
    end
end

% delete vertices of the outer approximation Tp of the upper image of (P)
% being very close to each other (for stability reasons) 
Tp=deletemultiples(Tp,options.close_vert_eps);

end

