function C = deletezeros(A,tol)
% delete zero columns in a matrix A
[m,n]=size(A);
B=zeros(m,n);
k=0;
for i=1:n
    if norm(A(:,i))>tol
        k=k+1;
        B(:,k)=A(:,i);
    end
end
C=B(:,1:k);
end

