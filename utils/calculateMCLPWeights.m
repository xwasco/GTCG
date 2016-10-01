function [weights,weightsTmp]=calculateMCLPWeights(As)

    S=[];
    for i=1:size(As,1)
        S(:,:,i)=As{i};
    end

    %This one combines using different weights and give Accuracy and the class

    n=size(S,1);
    weightsTmp=[];
    vv=size(S,3);
    U=zeros(vv,n);
    x=rand(n,1)+1000;
    x=x/sum(x); 

    for ih=1:size(S,3)
        U(ih,:)=x'*S(:,:,ih);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Call Bensolve using different options

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear options;
    options.lp_solver=3; %it is faster
    options.eps=1e-005; %stopping criteria
    options.alg3_eps=1e-007; 
    options.dual=1; %solve the dual problem, is simplerand since it a convex optimization the solution is the same for the primal and dual

    B=[ones(1,n);-ones(1,n);eye(n)]; 
    b=[1;-1;zeros(n,1)];
    %Pp=[U1;U2;U3;U4;U5;U6];
    Pp=U;
    k=size(Pp,1);
    [Shh,Sh,T,PP,PPh,DD]=bensolve(-Pp,B,b,[],[],[],options);
    weights=T(end-k+1:end,:);
    %weights

    jj=0;
    weightsTmp=[];
    for ni=1:size(weights,2)  
       ind=find(weights(:,ni)>=0.9999); %remove degenerate case (weight = 1 means selecting only 1 matrix)
       if isempty(ind)
           jj=jj+1;
           weightsTmp(:,jj)=weights(:,ni);
       end  
    end

    %sometime weights are very small but negative...fix this
    weightsTmp=abs(weightsTmp);
    
    if(jj~=0)
        weights=weightsTmp;
        weights=sum(weights,2)./sum(sum(weights));
    else
        %if no combination is non trivial set to equal weights
        weights=ones(size(weights,1),1)./size(weights,1);
        weightsTmp=weights;
    end     
    
end