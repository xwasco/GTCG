function [ A ] = zeropadding( A,index )
%ZEROPADDING Add a zero row and column in matrix A in correspondances of the
% index vector elements

index=sort(index);

for ix=1:numel(index)
    if index(ix)>1 && index(ix)<=size(A,1) 
        %add the row
        A=[A(1:index(ix)-1,:) ; zeros(1,size(A,2)) ; A(index(ix):end,:)];

        %add the column
        A=[A(:,1:index(ix)-1,:) , zeros(size(A,1),1) , A(:,index(ix):end)];
    elseif index(ix)>size(A,1)
        A=[A , zeros(size(A,1),1) ; zeros(1,size(A,1)+1)];
    else
        A=[zeros(1,size(A,2)+1) ; zeros(size(A,1),1) , A];
    end
    
    index=index+1; %update the index so
end

end

