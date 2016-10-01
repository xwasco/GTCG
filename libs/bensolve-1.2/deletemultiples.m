function C = deletemultiples(A , tol)
% delete multiple columns in a matrix
[m,n]=size(A);
B=zeros(m,n);
k=0;
for i=1:n
    if ~vectorismember(A(:,i),B(:,1:k),tol)
        k=k+1;
        B(:,k)=A(:,i);
    end
end
C=B(:,1:k);
end

