function B = normalize(A)
%NORMALIZE returns normalized column vectors of a matrix (w.r.t. 2-norm)
B=A./(ones(size(A,1),1)*sqrt(sum(A.^2,1)));
end

