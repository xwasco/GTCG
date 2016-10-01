function tf = vectorismember(x, S, tol)
% VECTORISMEMBER tests if a vector x is a column vector of the matrix S
% with tolerance tol
tf=0; 
for i=1:size(S,2)
    if norm(x-S(:,i)) < tol
        tf=1;
        break;
    end
end
end

