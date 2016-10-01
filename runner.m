% Sebastiano Vascon      Version 1.01
% Copyright 2016 Sebastiano Vascon.  [sebastiano.vascon-at-iit.it]
% Please email me if you have any questions.
%
% Please cite one of these works
% [1] S. Vascon, E. Zemene , M. Cristani, H. Hung, M.Pelillo and V. Murino 
% Detecting conversational groups in images and sequences: A robust game-theoretic approach
% Computer Vision and Image Understanding (CVIU), 2016
%
% [2] S. Vascon, E. Zemene , M. Cristani, H. Hung, M.Pelillo and V. Murino
% A Game-Theoretic Probabilistic Approach for Detecting Conversational Groups
% ACCV 2014

% -------------------------------------------- %
%                Detect Groups in Single Frame %
% -------------------------------------------- %

clc
clear
warning off;

addpath(genpath('utils'));
addpath(genpath('libs'));

%% ALGORITHM PARAMETERS

% dataset directory
datasetDir='data/Synth';
run([datasetDir '/dsetParameter.m']); %load the dataset parameters
param.evalMethod='card';                    %'card' require that 2/3 of individals are correctly matched in a group
                                            %'all'  require that 3/3 of individals are correctly matched in a group (more stricter evaluation)

%multi/single frame
param.numFrames=1;                          %number of frames to analyze (>1 imply multiframe analysis)

%parameter for the 2D histogram
param.hist.n_x=20;                          %number of rows for the frustum descriptor
param.hist.n_y=20;                          %number of columns for the frustum descriptor

%displaying options
param.show.weights=0;                       %show the weight used to condense the similarity matrices
param.show.groups=0;                        %show a figure with the current frame, the decetion and the groundtruth
param.show.frustum=0;                       %show the frustum 
param.show.results=1;                       %display the precision/recall/F1-score values

%weight calculation parameters
param.weight.mode='MOLP';                   %the multiframe mode is activated only if param.numFrames>1. Set to:
                                            %'MOLP' (MultiObjectiveLinearProgramming)
                                            %'EQUAL' use equal weights for the frames
                                            %'MAXENTROPY' pick the
                                            %combination that maximize the entropy of the weight

%set the frustum modality
frustumMode='CVIU';                         %'CVIU' use the CVIU model (cite [1])
                                            %'ACCV' use the ACCV model (cite [2])

seqDir=''; %if a sub-sequence exists write the folder name here
datasetDir=[datasetDir seqDir];

load([datasetDir '/features.mat'],'features','timestamp');
load([datasetDir '/groundtruth.mat'],'GTgroups','GTtimestamp');

[~,indFeat] = intersect(timestamp,int64(GTtimestamp));
timestamp = timestamp(indFeat);
features  = features(indFeat);

precisions=[]; 
recalls=[];
TPs=[];
FPs=[];
FNs=[];

detections=[];

for f=1:numel(features)
    if ~isempty(features{f})
        
        feat=features(f:f+param.numFrames-1);                   %copy the frames
        
        fprintf(['******* Frames ' num2str(f:f+param.numFrames-1) ' *******\n']);
        
        [groups, frustums,weights]=detectGroups(feat,param);    %detect groups
        
        if param.show.weights>0 
            %display the weights
            figure(param.show.weights);
            bar(weights);
            title(['Weights used in Eq 8 of ACCV']);
        end
        
        detections=[detections ; {groups}];
        
        if param.show.groups>0
            fr=f+param.numFrames-1; %frame used as reference for the evaluation
            showGroups(fr,groups,GTgroups,param);
        end
        
        [p,r,tp,fp,fn] = evalgroups(groups,GTgroups(:,f+param.numFrames-1),param.evalMethod);
        
        %add the evaluation results to the queue
        precisions=[precisions ; p'];
        recalls=[recalls ; r'];
        TPs=[TPs ; tp'];
        FPs=[FPs ; fp'];
        FNs=[FNs ; fn'];
        
        showResults(precisions,recalls);

    end
end

