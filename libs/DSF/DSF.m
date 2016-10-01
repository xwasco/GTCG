% Sebastiano Vascon      Version 1.00
% Copyright 2014 Sebastiano Vascon.  [sebastiano.vascon-at-iit.it]
% Please email me if you have questions.
%
% Please cite this work
% [1] S. Vascon, E. Zemene , M. Cristani, H. Hung, M.Pelillo and V. Murino
% A Game-Theoretic Probabilistic Approach for Detecting Conversational Groups
% ACCV 2014

classdef DSF
    %Dominant Set
    
    properties
        maxIters=20000; %Stop criteria: max number of iterations to converge to a cluster
        eps=1e-10;      %Stop criteria: distance between two succeded step
        noise=0;        %noise added to the initial vector
        A=[];           %payoff matrix
        dyn='RD';       %dynamic system to apply RD=replicator dynamics
        X=[];           %the dataset
        RD_exp_k=1;     %exponential parameter for the exponential replicator dynamics
        clusters=[];    %the resulting clusters set
        stat=[];        %some time statistics
        V=[];           
        startSize=0;    %population size
        globalContext=0;%stopping criterion proposed by Hung & Krose in "Detecting F-formations as Dominant Sets" ICMI 2011
        
        alternativeWeight=0; %alternative way of computing weight based on eq (6) of http://www.dais.unive.it/~pelillo/papers/PAMI%202007.pdf
        
        simFunct;       %similarity function. e.g. simFunct=@(x,y,f_param)exp(-sqrt(sum(x-y).^2)/f_param); to generate a classical gaussian kernel on the data
        
        W=[];
    end
    
    methods(Access=public)
       
        function obj=DSF(dyn,maxIters,eps,noise)
            % DSF  Object constructor.
            %   DYN Dynamic system to apply 'RD'=replicator dynamics,
            %   'II'=Infection Immunization ('II')
            %   MAXITERS Stop criteria: max number of iterations (5000)
            %   EPS Stop criteria: distance between two succeded step (1e-7)
            %   NOISE Noise to be added on the initialization phase.
            %   bibfn the filename in which the bibliography of the methods
            %   will be stored
            
            obj.dyn=dyn;
            obj.eps=eps;
            obj.maxIters=maxIters;
            obj.noise=noise;
        end
        
        function obj=setDataset(obj,X)
            obj.X=X;
        end
        
        function obj=setSimilarityMatrix(obj,A)
            obj.A=A;
            obj.startSize=size(A,1);
        end
        
        function obj=clusterize(obj,th)
            % Extract the clusters using the Dominant Set
            %   -th     represent the threshold to be applied at convergence to
            %           retrieve the "survived" set of strategies. Th represent the
            %           percentage over the best strategy (e.g. x=[.7 .2 .05 .05]
            %           th=.7 survided strategy are x_i>th*max(x)=.7*.7 => [x1]
            %           th=.1 survided strategy are x_i>th*max(x)=.1*.7 => [x1,x2]
            
            x=ones(size(obj.A,1),1)./size(obj.A,1)+rand(1)*obj.noise;
            x=x./sum(x);            
            obj=clusterizeBoosted(obj,th,x);
            
            %fprintf('Generated vector x\n');
            %obj=getClusters(obj,th,x);
        end
        
        function obj=clusterizeBoosted(obj,th,x)
            % Extract the clusters using the Dominant Set 
            %   -th     represent the threshold to be applied at convergence to
            %           retrieve the "survived" set of strategies. Th represent the
            %           percentage over the best strategy (e.g. x=[.7 .2 .05 .05]
            %           th=.7 survided strategy are x_i>th*max(x)=.7*.7 => [x1]
            %           th=.1 survided strategy are x_i>th*max(x)=.1*.7 => [x1,x2]
            %   -x      Initialize the starting solution
            
            if size(obj.A,1)>=4000 && (strcmp(obj.dyn,'RD')==1 || strcmp(obj.dyn,'RDEXP')==1 )
                fprintf('Since the size of the matrix is high (>4000) is suggested to use the Infection Imunization dynamics obj.dyn=''II''\n');
            end
            
            tA=obj.A;
            obj=getClusters(obj,th,x);
            obj.A=tA;
        end
        
        function [obj]=getDataCluster(obj)
        % Return a matrix of cells in which the first column is the cluster
        % and the second is a cell containing the original dataset entry.
            for i=1:numel(obj.clusters)
                obj.clusters(i).elements=obj.X(obj.clusters(i).index,:);
            end
        end
        
        function [idx]=getClusterDataAssociation(obj)
            % Return a column vector in which each cell represent the
            % cluster associated to the i-th element. i.e. idx(12)=cluster
            % to whom the 12-th element belong
            idx=zeros(size(obj.A,1),1);
            for i=1:numel(obj.clusters)
                for j=1:numel(obj.clusters(i).index)
                    idx(obj.clusters(i).index(j))=i;
                end
            end
        end
         
        function n=getNumberOfClusters(obj)
            % Get the number of clusters
            n=numel(obj.clusters);
        end
        
        function CA=getClustersAssignment(obj)
            CA=zeros(size(obj.A,1),size(obj.A,1));
            for c=1:numel(obj.clusters)
                for i=1:numel(obj.clusters(c).index)
                    for j=i+1:numel(obj.clusters(c).index)
                        CA(obj.clusters(c).index(i),obj.clusters(c).index(j))=1;
                        CA(obj.clusters(c).index(j),obj.clusters(c).index(i))=1;
                    end
                end 
            end
        end
        
        %calculates the average weighted degree of node wrt set
        function w=CalcAwdeg(obj,nodeSet,nodeId)
            w=0;
            
            if numel(nodeSet)==0
                fprintf('nodeSet size is 0!\n');
            else
                for i=1:numel(nodeSet)
                    w=w+obj.W(nodeId,nodeSet(i));
                end
                w=w/numel(nodeSet); 
            end
        end

        function w=CalcRelSim(obj,nodeSet,i,j)
            w=obj.W(i,j)-obj.CalcAwdeg(nodeSet,i);
        end

        function w=CalcNodeWeight(obj,currset,i,w)
            if obj.alternativeWeight==0
                if numel(currset)==1
                    w=1;
                elseif numel(currset)==2 
                    w=obj.W(currset(1),currset(2));
                else
                    newset=currset;
                    id=find((newset==i)==1);
                    newset(id)=[]; %S\{i}
                    weights=zeros(numel(newset),1);
                    for j=1:numel(newset)
                        temprelsim=obj.CalcRelSim(newset,newset(j),i);
                        tempweight=obj.CalcNodeWeight(newset,newset(j),w);
                        weights(j)=weights(j)+temprelsim*tempweight;
                    end
                    w=sum(weights);
                    newset=[];
                end
            else
                if numel(currset)==1
                    w=1;
                else 
                    %choose an h in S\{i}
                    w=0;
                    S=currset;
                    %fprintf(['\nS={' num2str(S) '} , i=' num2str(i) ' => ']);
                    S(find(currset==i))=[];  
                    %fprintf(['S\\{' num2str(i) '} = {' num2str(S) '}\n']);
                    
                    for c=1:numel(S)
                        h=S(1);
                        j=S(c);                        
                        %fprintf(['(A(' num2str(i) ',' num2str(j) ')-A(' num2str(h) ',' num2str(j) ')*']);
                        %fprintf([num2str(S) ' - ' num2str(i) ' : ' num2str(h) '\n']);
                        %if h~=j
                            w=w+(obj.W(i,j)-obj.W(h,j))*obj.CalcNodeWeight(S,j,w);
                        %end
                    end                    
                end
            end
        end
        
        function [w]=iterativeCalcNodeWeight(S,i)
            w=0;
            currSet=S;
            currSet(find(currSet==i))=[];
            while numel(currSet)>1
                t=randperm(numel(S));
                h=S(t(1));
                w=w+(obj.W);
            end
            w=w+1;
        end
        
    end
    
    methods (Access=private)

        function [W]=getSimilarityMatrix(obj,X,funct,f_param,sim)
            % GENERATESIMILARITYMATRIX Generate the similarity matrix used
            % to cluster the dataset
            % INPUT:
            %   X    The row wise dataset (X(i,j)= the j-th propoerty
            %           of the i-th item
            %   funct   Handle for the function that calculate the pairwise
            %           similarity value. i.e.
            %           funct=@(x,y,f_param)exp(-sqrt(sum(x-y).^2)/f_param)
            %   sim     sim=1 the similarity matrix is simmetric (half
            %           computation time), sim=0 the matrix is asymmetric
            
            W=zeros(size(X,1),size(X,1));
            
            if sim==1
                for i=1:size(W,1)
                    for j=i+1:size(W,1)
                        W(i,j)=funct(X(i,:),X(j,:),f_param);
                    end
                    fprintf([num2str(i) '\n']);
                end
                W=W+W';
            else
                for i=1:size(W,1)
                    for j=1:size(W,1)
                        W(i,j)=feval(funct, X(i,:),X(j,:),f_param);
                        %A(j,i)=A(i,j);
                    end
                end
                W=W.*not(eye(size(W,1),size(W,1))); %make the diagonal equal to zero
            end
        end
        
        function [W]=getSimilarityMatrixWholeSet(obj,X,funct,f_param)
            % GENERATESIMILARITYMATRIX Generate the similarity matrix used
            % to cluster the dataset. In this case the funct accept as
            % input the whole dataset and a param (which could be a struct
            % with additional parameter).
            % INPUT:
            %   X    The row wise dataset (X(i,j)= the j-th propoerty
            %           of the i-th item
            %   funct   Handle for the function that calculate the pairwise
            %           similarity value. i.e.
            %           funct=@(x,f_param)exp(-sqrt(sum(x-y).^2)/f_param)
            %   sim     sim=1 the similarity matrix is simmetric (half
            %           computation time), sim=0 the matrix is asymmetric
            
            W=zeros(size(X,1),size(X,1));
            for i=1:size(W,1)
                    for j=i+1:size(W,1)
                        W(i,j)=feval(funct, X(i,:),X(j,:),f_param);
                        W(j,i)=W(i,j);
                    end
            end
            [~ ,W]=funct(X,f_param);
            W=W.*not(eye(size(W,1),size(W,1))); %make the diagonal equal to zero
        end
        
        function obj=getClusters(obj,th,x)
            % GETCLUSTERS  Extract the cluster 
            %   TH represent the threshold to be applied at convergence to
            %   retrieve the "survived" set of strategies.
            %   In property "clusters" an array of clusters will be created
            %   X represent the initialization vector
            %   tA is the original similarity matrix
            
            x=x./sum(x); %guarantee that is normalized
            
            obj.clusters=[];
            dict=1:size(obj.A,1);
            
            popSize=size(obj.A,1);
            
            obj.stat.maxSize=-inf;
            obj.stat.minSize=inf;
            obj.stat.elapsedTime = cputime;
            obj.V=zeros(popSize,1);
            obj.W=obj.A;
            
            StillGC=1;
            while(size(x,1)>0 && StillGC) %till there is something to extract
                fprintf(['Extraction of cluster ' num2str(numel(obj.clusters)+1) ' started... ']);
                
                c=getCluster(obj,x,th);
                
                %% MOD SEBA cohesiveness = 0 
                if c.cohesiveness<eps || isnan(c.cohesiveness)
                    c.cohesiveness=0; %sum(sum(x.*((obj.A.*not(eye(size(obj.A))))*x)));
                    c.population=zeros(numel(x),1);
                    c.index=find(x>0);
                    c.index=c.index(1);
                    c.populationThresholded=zeros(numel(x),1);
                    c.threshold=0;
                    c.stat.iteration=0;
                    c.stat.error=0;
                    c.centroid.value=0;
                    c.centroid.id=1;
                end
                
                %if global context is active use it to stop the extraction
                %of the clusters
                if(obj.globalContext==1)
                    % check the weight of each node not in DS that if added
                    % should be negative 
                    
                    w=0;
                    
                    ids=dict(c.index); 
                    
                    
                    %ids=c.index';
                     
                    for i=1:size(obj.W,1)
                        %add the i-th node into the set and evaluate its
                        %weight
                        
                        if isempty(find(ids==i))
                            w=0;
                            % if node is not into the DS evaluate the 
                            % weight to see the consistency of the DS of being 
                            % dominant
                            %w=obj.CalcNodeWeight([ids i],i,w);
                            w=CalcWeight(obj.W,[ids i],i); %using MEX FILE
                            %fprintf(['\n' 'w([' num2str(ids) '],' num2str(i) ')=' num2str(w)]);
                            if w>=0
                                StillGC=0;
                                break;
                            end
                        end
                    end
                end
                
                if obj.globalContext==0 || (obj.globalContext==1 && StillGC==1)
                    
                    %update some statistics
                    if obj.stat.minSize>size(c.index,1)
                        obj.stat.minSize=size(c.index,1);
                    end

                    if obj.stat.maxSize<size(c.index,1)
                        obj.stat.maxSize=size(c.index,1);
                    end

                    %erase the rows and column so the problem get smaller at
                    %each iteration
                    obj.A(c.index,:)=[];
                    obj.A(:,c.index)=[];

                    id=c.index;

                    obj.V(dict)=obj.V(dict)+c.population;
                    population=zeros(obj.startSize,1);
                    population(dict(c.index))=c.population(c.index);

                    c.populationThresholded=population;

                    c.index=dict(c.index);
                    c.centroid.id=dict(c.centroid.id);

                    obj.clusters=[obj.clusters c];

                    dict(id)=[];

                    x=ones(size(obj.A,1),1)./size(obj.A,1)+rand(1)*obj.noise;
                    x=x./sum(x);
                end
                fprintf(['done (#iter=' num2str(c.stat.iteration) ' , #err=' num2str(c.stat.error) ') !\n']);                
                
            end
            obj.V=obj.V./sum(obj.V);
            
            obj.stat.elapsedTime = cputime-obj.stat.elapsedTime;
        
        end
        
        function clust=getCluster(obj,x,th)
            dist=inf;
            iter=0;
            
            clust.index=[];
            clust.centroid.id=0;
            clust.centroid.value=0;
            clust.cohesiveness=0;
            
            if (strcmp(obj.dyn,'RD')==1)
                while (obj.eps<dist && iter<obj.maxIters)
                    old_x=x;
                    x=x.*(obj.A*x);
                    x=x./sum(x);
                    dist= pdist([x,old_x]');
                    iter=iter+1;
                end
            else
                [iter , dist,  pop] = iidyn(obj.A, x', obj.eps, obj.maxIters);
                x=pop';
            end
            
            %clust.cohesiveness=sum(sum(x.*(obj.A*x)));
            
            clust.cohesiveness=sum(sum(x.*((obj.A.*not(eye(size(obj.A))))*x)));

            
            t=max(x)*th;
            clust.population=x;
            clust.index=find(x.*(x>(t)));
            clust.populationThresholded=x.*(x>(t));
            clust.threshold=t;
            clust.stat.iteration=iter;
            clust.stat.error=dist;
            [clust.centroid.value , clust.centroid.id]=max(x);
        end
        
        function x=getComposedPopulation(obj)
            tot=0; %population size
            for i=1:obj.getNumberOfClusters()
                tot=tot+size(obj.clusters(i).index,2);
            end
            
            x=zeros(1,tot);
            for i=1:tot
                for j=1:obj.getNumberOfClusters()
                    x(i)=x(i)+obj.clusters(j).population(i);
                end
            end
        end

    end

end

