function [ ptsFM , orj,len] = frustumModel( pos,orj,length,aperture,samples,mode)
% FRUSTUMMODEL generate samples distributed on the field of view
% 
% INPUT:
% ======
%   pos         :=  position of a person [x,y]
%
%   orj         :=  head orientation between 0 and 2*pi
%
%   length      :=  length of the field of view
%
%   aperture    :=  aperture of the field of view in degree, for human 
%                   is 160 (default 160).
%
%   samples     :=  number of samples to generate (default 2000)
%
%   mode        :=  'CVIU' or 'ACCV' (default) 'beta'
% 
% OUTPUT:
% =======
%	ptsFM       :=  the samples of the frustum
%
%
% EXAMPLE:
% =======
%   Generate a frustum of 5000 points for a person in x=0 , y=0 with head
%   orientation of pi/2 (90 degree) with a length of view equal to 60 and
%   an aperture of 160 degree.
%
%   figure; fpts=frustumModel([0,0],pi/2,60,160,5000); scatter(fpts(:,1),fpts(:,2));
%   figure; fpts=frustumModel([0,0],pi/2,60,160,5000); surf(hist2D(fpts(:,1),fpts(:,2),20,20)); shading interp;
%
% -------------------------------------------------------------------------
% Sebastiano Vascon      Version 1.00
% Copyright 2014 Sebastiano Vascon.  [sebastiano.vascon-at-iit.it]
% Please email me if you have questions.
%
% Please cite this work if mode ACCV is used
% A Game-Theoretic Probabilistic Approach for Detecting Conversational Groups
% S. Vascon, E. Zemene , M. Cristani, H. Hung, M.Pelillo and V. Murino
% ACCV 2014
%
% Please cite this work if mode CVIU is used
% Detecting conversational groups in images and sequences: A robust game-theoretic approach
% S. Vascon, E. Zemene , M. Cristani, H. Hung, M.Pelillo and V. Murino
% CVIU 2016 doi:10.1016/j.cviu.2015.09.012

if nargin<6
    mode='ACCV';
end

if nargin<5
    samples=2000;
end

if nargin<4
    aperture=160;
end

pts=[];
        
ptsFM=[];
        
aperture=aperture/2;

if length>0
    
    if strcmp(mode,'CVIU')==1
        %length=length*3;
        %generate the random orientations
        sigma=((aperture/360*(2*pi)))/3;
        orj=normrnd(orj,sigma,samples,1);
        len=betarnd(0.8,1.1,samples,1).*length;
        len=sort(len,'ascend');
        ptsFM=[pos(1)+cos(orj).*len , pos(2)+sin(orj).*len];
    elseif strcmp(mode,'beta')==1
        %length=length*3;
        %generate the random orientations
        sigma=((aperture/360*(2*pi)))/3;
        orj=normrnd(orj,sigma,samples,1);
        len=betarnd(1.2,1.1,samples,1).*length;
        len=sort(len,'ascend');
        ptsFM=[pos(1)+cos(orj).*len , pos(2)+sin(orj).*len];
    elseif strcmp(mode,'CVIUREB')==1
        %generate the random orientations
        sigma=sqrt(((aperture/360*(2*pi)))/6);
        orjs=normrnd(orj,sigma,samples,1);
        len=betarnd(0.8,1.1,samples,1).*length;
        len=sort(len,'ascend');
        ptsFM=[pos(1)+cos(orjs).*len , pos(2)+sin(orjs).*len];        
    elseif strcmp(mode,'CVIU_Savarese')==1
        %length=length*3;
        %generate the random orientations
        sigma=((aperture/360*(2*pi)))/3;
        orj=normrnd(orj,sigma,samples,1);
        
        len=betarnd(0.8,1.1,samples,1).*length;
        ptsFM=[pos(1)+cos(orj).*len , pos(2)-sin(orj).*len];        
    elseif strcmp(mode,'CVIUGamma')==1
        %length=length*3;
        %generate the random orientations
        sigma=((aperture/360*(2*pi)))/3;
        orj=normrnd(orj,sigma,samples,1);
        len=gamrnd(1,4,samples,1).*length;
        ptsFM=[pos(1)+cos(orj).*len , pos(2)+sin(orj).*len]; 
    else
        orjX=cos(orj)*length;
        orjY=sin(orj)*length;

        %convert orj in deg
        orj=360*orj/(2*pi);

        % remove the point not in the frustum
        while size(ptsFM,1)<samples
            pts=[];
            for i=1:size(pos,2)
                pts=[pts , normrnd(0,length,samples,1)];
            end
            if ~isempty(pts)
                for i=1:size(pts,1)

                    pt=pts(i,:);

                    dot=sum(pt.*[orjX orjY]);
                    L=sqrt(sum(pt.^2));

                    ang=acos(dot/(L*length))/(2*pi)*360;% compute the angle between the point and the orientation

                    if ang<=(aperture/2) && size(ptsFM,1)<samples
                        %the point is inside the cone
                        ptsFM=[ptsFM ; pts(i,:)+pos];
                    end
                end
            end
        end
    end
else
    ptsFM=repmat(pos,samples,1);
end
end
