function [groups,px,weights,ALLWeights]=detectGroups(frames,param)
% DETECTGROUPS Detect conversational groups given a set of frames each one 
% containing the positions, the head orientation and the IDs of the 
% individuals in the scene.
%
% The last frame is the reference frame because the method used the past
% information to smooth the actual values. For this reason the further 
% evaluation should be made on it.
% 
% INPUT:
% ======
%  frames     := a 1xm vector of cells containing on each cell a matrix of
%                persons. Each matrix is of size nx4 having on each row: 
%               | Person ID | Pos X | Pos Y | Orientation (radiants) |
%
%  param      := a struct containing the following parameters. If not
%                specified a default set of parameters will be used.
%
%                Frustum parameters:
%                   - param.frustum.length=20;      lenght of the frustum
%                   - param.frustum.aperture=160;   aperture in degree
%                   - param.frustum.samples=2000;   number of samples
%                   - param.frustum.mode='CVIU';    CVIU or ACCV
%
%                Histogram paramters
%                   - param.hist.n_x=20;            number of bin per rows
%                   - param.hist.n_y=20;            number of bin per cols    
%                Affinity matrix parameter
%                   - param.sigma=0.2;              normalizing factor 
%                   - param.checkFacing=0;          check if the persons are facing
%                   - param.method='JS';            type of distance (JS or
%                                                   KL)
%                Weights for multi-frames
%                   - param.weight.mode;            if not set weights are found by using MLCP 
%                                                   'EQUAL' use equal weights for the frames
%                                                   'MAXENTROPY' use the
%                                                   set of weights that
%                                                   maximize the entropy
%                                                   (the set close to equal
%                                                   weight)
%
% OUTPUT:
% =======
%  groups     := vector of cells in which each cell represents a group 
%                containing the Persons ID 
%
%  frustums   := a list of cells containing the generated frustum per each
%                persons and the corresponding 2D histogram
%
%  weights    := the average of the weights that the MCLP optimizer has 
%                assigned to each frame
%
%  ALLWeights := the base set of weights generated
%
%
% EXAMPLE:
% =======
%  frames=[];
%  frame_t1=[1 , 100, 100 , pi; 2 , 120, 120 , pi/2; 3 , 90, 90 , pi/4]; 
%  frames=[frames , {frame_t1}];
%  frame_t2=[1 , 110, 110 , pi; 2 , 125, 125 , pi/2; 3 , 100, 100 , pi/4]; 
%  frames=[frames , {frame_t2}];
%  frame_t3=[1 , 105, 110 , pi; 2 , 130, 125 , pi/2; 3 , 105, 95 , pi/4]; 
%  frames=[frames , {frame_t3}];
%  [groups, frustums,weights]= detectGroups(frames);
% 
%
% -------------------------------------------------------------------------
% Sebastiano Vascon      Version 1.01
% Copyright 2016 Sebastiano Vascon.  [sebastiano.vascon-at-iit.it]
% Please email me if you have any questions.
%
% Please cite one (or both) of these works
% [1] S. Vascon, E. Zemene , M. Cristani, H. Hung, M.Pelillo and V. Murino 
% Detecting conversational groups in images and sequences: A robust game-theoretic approach
% Computer Vision and Image Understanding (CVIU), 2016
%
% [2] S. Vascon, E. Zemene , M. Cristani, H. Hung, M.Pelillo and V. Murino
% A Game-Theoretic Probabilistic Approach for Detecting Conversational Groups
% ACCV 2014
% -------------------------------------------------------------------------

if nargin<2
    %frustum
    param.frustum.length=50;
    param.frustum.aperture=160;
    param.frustum.samples=2000;

    % histogram binning
    param.hist.n_x=20;
    param.hist.n_y=20;

    % affinity matrix parameter
    param.sigma=0.4;
    param.method='JS';
    param.checkFacing=1;
 
    %head orientation quantization
    param.HO_quantization=0;
end

if ~isfield(param.frustum,'mode') 
    param.frustum.mode='CVIU';
end

%find the set of IDs on all the frames and zero pad for the ones that
%are not in the list.
IDs=[];
for f=1:numel(frames)
    frame=frames{f};
    [IDs]=[IDs ; frame(:,1)];
    frame=sortrows(frame,1); %sort the ID the ensure consistency when zero 
                             %pad is added to the final matrix
    frames{f}=frame;
end
IDs=unique(IDs);    %find the unique ID
IDs=sort(IDs);      %sort the IDs for later use


if isfield(param, 'FillMissDetection') &&  param.FillMissDetection>0
    %in case of miss detection it fill the person using the mean of the
    %other frames
    mds={};
    for f=1:numel(frames)
        %find the miss detection
        mds{f}=find(ismember(IDs,frames{f}(:,1))==0);
    end
    
    for f=1:numel(mds)
        
        %search for the frames that do not have the missdetection of frame
        %f
        
        md=IDs(mds{f});
        if ~isempty(md)
            for j=1:numel(md)
                nmds={};
                for m=1:numel(mds)
                    if f~=m
                       nmds{m}=find(ismember(frames{m}(:,1),md(j))==1);
                    else
                       nmds{m}=[];
                    end
                end
                
                %in nmds there are the pointer to the data in the other frames,
                %now we make an average and substitute the miss detection
                tData=zeros(1,size(frames{f},2));
                counter=0;
                for i=1:numel(nmds)
                    if ~isempty(nmds{i})
                        tData=tData+frames{i}(nmds{i},:);
                        counter=counter+1;
                    end
                end
                tData=tData./counter;
                frames{f}=[frames{f} ; tData];
                frame=sortrows(frames{f},1);
                frames{f}=frame;
            end
        end
    end 
    %{
    for f=1:numel(frames)
        tFrame=[];
        for p=1:numel(IDs)
            for i=1:size(frames{f},1)
                if frames{f}(i,1)==IDs(p)
                    tFrame=[tFrame ; frames{f}(i,:)]; 
                else
                    %copy the average informations from the other frames
                    
                    tFrame=[tFrame ; t];
                end
            end
        end
        frames{f}=tFrame;
        [IDs]=[IDs ; frame(:,1)];
        frame=sortrows(frame,1); %sort the ID the ensure consistency when zero 
                                 %pad is added to the final matrix
        frames{f}=frame;
    end
    
    r={};
    for i=1:numel(As)
        r{i}=find(sum(As{i},2)==0);
    end
    
    %evaluate average filling for each empty row from the other matrices
    for i=1:numel(As)
        if ~isempty(r{i}) %if there are some row to fill in the matrix As{i}
            for j=1:numel(r{i}) %for each empty row
                rF=zeros(1,size(As{1},1));
                cF=zeros(size(As{1},1),1);
                counter=0;
                rId=r{i}(j);
                for k=1:numel(As) %for all matrices except the i-th and the one that has the same empty row
                    if k~=i && ~ismember(rId,r{k})
                        rF=rF+As{k}(rId,:);
                        cF=cF+As{k}(:,rId);
                        counter=counter+1;
                    end
                end
                if counter==0
                    fprintf('Err');
                end
                rF=rF./counter;
                cF=cF./counter;
                As{i}(rId,:)=rF;
                As{i}(:,rId)=cF;
            end
        end
    end
    %}
end


%generate a similarity matrix for each frame leaving the space for the
%frame in which persons are not present.
As=[];
for f=1:numel(frames)

    persons=frames{f};
    
    % generate a frustum for each person
    minX=inf;
    minY=inf;
    maxX=-inf;
    maxY=-inf;

    px={};
    for i=1:size(persons,1)
        
        if isfield(param,'HO_quantization') && param.HO_quantization>0
            px{i}.frustum=frustumModel([persons(i,2),persons(i,3)],angleQuantization(persons(i,4),[0,2*pi],param.HO_quantization),param.frustum.length,param.frustum.aperture,param.frustum.samples,param.frustum.mode);
        else
            px{i}.frustum=frustumModel([persons(i,2),persons(i,3)],persons(i,4),param.frustum.length,param.frustum.aperture,param.frustum.samples,param.frustum.mode);
        end
        
        %evaluate space boundaries for the histogram
        if minX>min(px{i}.frustum(:,1)) 
            minX=min(px{i}.frustum(:,1));
        end
        if minY>min(px{i}.frustum(:,2)) 
            minY=min(px{i}.frustum(:,2));
        end
        if maxX<max(px{i}.frustum(:,1)) 
            maxX=max(px{i}.frustum(:,1));
        end
        if maxY<max(px{i}.frustum(:,2)) 
            maxY=max(px{i}.frustum(:,2));
        end

    end

    %generate 2D histogram for each person and on the entire space
    for i=1:size(persons,1)
        px{i}.hist2D=hist2D(px{i}.frustum(:,1),px{i}.frustum(:,2),param.hist.n_x,param.hist.n_y,[minX,maxX],[minY,maxY]); %get the 2D histogram
        px{i}.hist=reshape(px{i}.hist2D,1,param.hist.n_x*param.hist.n_y); %create a row vector of the histogram
    end
    
    %evaluate pairwise affinity matrix    
    if numel(px)>1
        A=zeros(numel(px),numel(px));
        if strcmp(param.method,'JS')==1
            for i=1:numel(px)
                for j=i+1:numel(px)
                    A(i,j)=JSDiv(px{i}.hist,px{j}.hist);
                    A(j,i)=A(i,j);
                end
            end
        elseif strcmp(param.method,'KL')==1
            for i=1:numel(px)
                for j=1:numel(px)
                    if i~=j
                        A(i,j)=KLDiv(px{i}.hist,px{j}.hist);
                    end
                end
            end    
        else
            fprintf('Unrecognized distance function (use JS or KL)\n');
        end
        A=exp(-A./param.sigma).*not(eye(size(A)));
    else
        A=0;
    end
    
    %zero padding the matrix so all the matrices have the same size
    if size(A,1)<numel(IDs)
        %find the row and cols to pad
        ind=find(ismember(IDs,persons(:,1))==0);
        
        A=zeropadding(A,ind);
    end
    As=[As ; {A}];
end

EnsembleWeights=[];

if numel(As)>1
    %if the number of frames is greater than one we are working in a
    %multiframe setting and thus we have to compute the weight to mix 
    %the affinity matrices
    
    %find the weights using MCLP
    weights=[];
    if isfield(param,'weight') && isfield(param.weight,'mode')
        if strcmp(param.weight.mode,'EQUAL')
            weights=ones(numel(As),1)./numel(As);
        elseif strcmp(param.weight.mode,'MAXENTROPY')
            [~,weights] =calculateMCLPWeights(As);
            
            ALLWeights=weights;
            
            %calculate entropy for each set of weights
            Hs=zeros(1,size(weights,2));
            for wi=1:size(weights,2)
                Hs(wi)=weights_entropy(weights(:,wi));
            end
            
            [v,indH]=max(Hs);
            
            weights=weights(:,indH);
            
            fprintf(['Maximal entropy ' num2str(v) ' correspondes to weights: ' num2str(weights') '\n']);
        elseif strcmp(param.weight.mode,'ENSEMBLE')
            [weights,EnsembleWeights] =calculateMCLPWeights(As);  
            ALLWeights=EnsembleWeights;
            
        else
            [weights,EnsembleWeights] =calculateMCLPWeights(As); 
            ALLWeights=EnsembleWeights;
        end
    else
        [weights,ALLWeights]=calculateMCLPWeights(As); 
    end
else
    A=As{1};
    weights=1;
    ALLWeights=1;
    EnsembleWeights=ALLWeights;
end

if strcmp(param.weight.mode,'MOLP_DS')
    if size(ALLWeights,2)>2
        WA=pdist(ALLWeights');
        WA=squareform(exp(-WA));
        d=DSF('RD',1000,1e-5,1e-3);
        d=d.setSimilarityMatrix(WA);
        d=d.clusterize(0.0001);
        %d=d.setDataset(ALLWeights');
        %d=d.getDataCluster();
        t=sum(repmat(d.clusters(1).population',size(ALLWeights,1),1).*ALLWeights,2);
        weights=t./sum(t);
        %weights=ALLWeights(:,d.clusters(1).centroid.id);
    else
        weights=ones(size(ALLWeights))./size(ALLWeights,1);
    end
end

if ~isfield(param,'weight') || ~isfield(param.weight,'mode') || ~strcmp(param.weight.mode,'ENSEMBLE')
    A=zeros(size(As{1}));
    for i=1:numel(As)
        A=A+weights(i).*As{i};
    end
        
    % extract groups using dominant set library
    d=DSF('RD',3000,1e-5,0);
    %d=DSF('RD',3000,1e-10,0); %%SEBA
    if ~isfield(param,'globalcontext')
        d.globalContext=1;
    else
        d.globalContext=param.globalcontext;
    end
    
    theta=0.00001;
    %theta=0.0000001;
    
    d=d.setSimilarityMatrix(A);
    d=d.clusterize(theta);
    d=d.setDataset(IDs);
    d=d.getDataCluster();

    % get groups
    groups=[];
    for i=1:size(d.clusters,2)
        if numel(d.clusters(i).elements)>1
            groups=[groups; {d.clusters(i).elements'}];
        end
    end
elseif strcmp(param.weight.mode,'ENSEMBLE')
    %for each set of weight compute the clusters and find an agreements
    %between the solutions
    groupsEnsemble=zeros(size(As{1}));
    
    for wi=1:size(EnsembleWeights,2)
        w=EnsembleWeights(:,wi);
        
        %create a pairwise matrix using the weighted average
        A=zeros(size(As{1}));
        for i=1:numel(As)
            A=A+w(i).*As{i};
        end
        
        % extract groups using dominant set library
        d=DSF('RD',3000,1e-5,0);
        if ~isfield(param,'globalcontext')
            d.globalContext=1;
        else
            d.globalContext=param.globalcontext;
        end
        theta=0.00001;
        %d.globalContext=0;
        d=d.setSimilarityMatrix(A);
        d=d.clusterize(theta);
        d=d.setDataset(IDs);
        d=d.getDataCluster();

        % get groups
        groups=[];
        for i=1:size(d.clusters,2)
            groups{i}=d.clusters(i).elements';
            for j=1:numel(d.clusters(i).index)
                for k=j+1:numel(d.clusters(i).index)
                    groupsEnsemble(d.clusters(i).index(j),d.clusters(i).index(k))=groupsEnsemble(d.clusters(i).index(j),d.clusters(i).index(k))+1;
                    groupsEnsemble(d.clusters(i).index(k),d.clusters(i).index(j))=groupsEnsemble(d.clusters(i).index(j),d.clusters(i).index(k));
                end
            end
        end
       
    end
    if size(EnsembleWeights,2)>1
        % searching for ensemble consistency
        fprintf(['Searching for the ensemble consistency among ' num2str(size(EnsembleWeights,2)) ' possible clustering results\n']);
        % extract groups using dominant set library
        d=DSF('RD',3000,1e-5,1e-4);
        d.globalContext=1;
        theta=0.00001;

        d=d.setSimilarityMatrix((groupsEnsemble./size(EnsembleWeights,2))); %.^size(EnsembleWeights,2));
        d=d.clusterize(theta);
        d=d.setDataset(IDs);
        d=d.getDataCluster();

        % get groups
        groups=[];
        for i=1:size(d.clusters,2)
            groups{i}=d.clusters(i).elements';
        end
    end
    fprintf(['Found ' num2str(numel(groups)) ' groups\n']); 
end
end