function [frames,fid,pid]=injectNoise(frames,err)
% INJECTNOISE Add gaussian noise to the head orientation of person in a
% set of frames
%
% INPUT:
% ======
%  frames     := a 1xm vector of cells containing on each cell a matrix of
%                persons. Each matrix is of size nx4 having on each row: 
%               | Person ID | Pos X | Pos Y | Orientation (radiants) |
%
%  err        := a struct containing the following parameters. If not
%                specified no noise will be added.

%                err.percFrame=0;     percentage of corrupted frames
%
%                err.percPersons=0;   percentage of person corrupted on 
%                                       each currupted frame
%
%                err.noiseAmount=0;   amount of noise to be added [0,2pi]
%
% OUTPUT:
% =======
%  frames     := a 1xm vector of cells containing on each cell a matrix of
%                persons altered with noise. 
%                Each matrix is of size nx4 having on each row: 
%                | Person ID | Pos X | Pos Y | Orientation (radiants) |
%
%  fid        := the id of the modified frames
%
%  pid        := the id of the modified persons for each modified frame
%
%
% -------------------------------------------------------------------------
% Sebastiano Vascon      Version 1.00
% Copyright 2014 Sebastiano Vascon.  [sebastiano.vascon-at-iit.it]
% Please email me if you have questions.
%
% Please cite this work
% [1] S. Vascon, E. Zemene , M. Cristani, H. Hung, M.Pelillo and V. Murino
% A Game-Theoretic Probabilistic Approach for Detecting Conversational Groups
% ACCV 2014
% -------------------------------------------------------------------------

fid=[]; pid=[];

if nargin<2
    err.percFrame=0;
    err.percPersons=0;
    err.noiseAmount=0;
end

if err.percFrame>0 && err.percPersons>0 && err.noiseAmount>0

    %select frame based on the percentage of corruption
    numFrames=ceil(numel(frames)*err.percFrame/100);
    fid=randperm(numel(frames)); fid=sort(fid(1:numFrames));
    
    %select persons for each each corrupted frames
    for f=1:numel(fid)
        pers=frames{fid(f)};    
        numPers=ceil(size(pers,1)*err.percPersons/100);
        cp=randperm(size(pers,1)); cp=sort(cp(1:numPers));
        pid=[pid ; {cp}];
        
        for p=1:numel(cp)
            pers(cp(p),4)=pers(cp(p),4)+normrnd(0,err.noiseAmount);
        end
        frames{fid(f)}=pers;
    end
end

end