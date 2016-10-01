function [ Tp ] = vert2( Td, Td_hat, d_hat , Y , c )
% VERT2 chooses a vertex enumeration method for the dual Benson algorithm
%
% REFERENCES: see bensolve.m 
% LICENSE: This file is part of BENSOLVE, see bensolve.m

global options

% delete vertices of the outer approximation Td of the upper image of (P)
% being very close to each other (for stability reasons)
Td=deletemultiples(Td,options.close_vert_eps);

if size(Td,1)==1 % the scalar case
    Tp=max(Td);
else 
    % Choose the variant of the vert2 subroutine 

    % based on convhulln and dual polyhedra
    if options.vert_enum=='A'
        Tp=vert2a(Td,Td_hat,Y,c);
        
    % based on convhulln and polyhedron-polytope transformation    
    elseif options.vert_enum=='B'
        Tp=vert2b(Td,Td_hat,d_hat,Y,c);
    
    % based on cdd and polyhedron-polytope transformation
    elseif options.vert_enum=='C'
        Tp=vert2c(Td,Td_hat,d_hat,Y,c);
    else
        error('Chosen variant of for vertex enumeration does not exist.');
    end
end

% delete vertices of the outer approximation Tp of the lower image of (D^*)
% being very close to each other (for stability reasons) 
Tp=deletemultiples(Tp,options.close_vert_dual_eps);

end