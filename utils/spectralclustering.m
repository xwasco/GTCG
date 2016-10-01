function [ d ] = spectralclustering( A,IDs)
%SPECTRALCLUSTERING Summary of this function goes here
%   Detailed explanation goes here

d.clusters=[];
idD=[];
if nargin<2
    idD=[1:size(A,1)]';
else
    idD=IDs;
end
    
if size(A,1)>1

    %A=A+double(A>0);
    %A=A+eye(size(A));
    %A=A+eps; %ones(size(A)); %eps;

    k_min=0;
    k_max=size(A,1);
    
    D=zeros(size(A,1));
    for i=1:1:size(A,1)
        D(i,i)=sum(A(i,:));
    end
    L= D-A;
    L=D^(-1/2)*L*D^(-1/2);

    [U,V]= eig(L);
    v=diag(V);
    [vs, is] = sort(v,'ascend');

    eigengaps = zeros(length(vs)-1,1);

    for i=1:1:length(eigengaps)
        if ((i<k_min) || (i> k_max))
            eigengaps(i)=-1;
        else
            eigengaps(i)=vs(i+1)-vs(i);
        end
    end

    [~ , k] = max(eigengaps);

    fprintf([' Found ' num2str(k) ' clusters\n']);

    % the number of eigenvectors that we will use is k
    % we now need to create a length(data(:,1)) x k matrix each column of which
    % are the k smallest eigenvector of the normalized Laplacian]
    u=[];
    for i=1:1:k 
        u = [u U(:,is(i))];
    end

    u=real(u);    
    % we re-normalize the matrix rows to unit norm 
    for i=1:1:length(u(:,1))
        u(i,:)=u(i,:)/sqrt(sum(u(i,:).*u(i,:)));
    end
    
    u=real(u);
    % we are using k-means.  You can use any other clustering algorithm 
    % we also use the parameter replicates in order to obtain the "steady
    % state" solution for k-means
    %if k<=u
    try
        [C, C_spec] = kmeans(u,k,'Replicates',200);
    catch ME
        fprintf(ME);
    end
    
    count=1;
    for i=1:1:k
        idx=find(C == i);
       % if numel(idx)>1 
            d.clusters(count).elements=idD(idx);
            count=count+1;
       % else
        %    fprintf('Group discarderd because composed only by one person\n');
       % end
    end
    
elseif size(A,1)==1
    d.clusters(1).elements=idD(1);
else
    
end

end

