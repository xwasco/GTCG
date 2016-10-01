function [ m,stepX,stepY ] = hist2D(x, y, nbins, mbins,xrange,yrange)
% HIST2D Generate a 2D histogram of a given set of points. The histogram is
% normalized so it adds up to 1.
%
% INPUT:
% ======
%   x        := a vector of mx1 position in the x directions
%
%   y        := a vector of mx1 position in the y directions
%
%   nbins    := number of bins (rows)
%
%   mbins    := number of bins (columns)
%
%   xrange   := [minX,maxX] a vector of two elements that define the range
%               in x direction
%
%   yrange   := [minY,maxY] a vector of two elements that define the range
%               in y direction
%
% OUTPUT:
% =======
%  m         := a nbins x mbins matrix of double containing the 2D
%               histogram
%
%  stepX     := the step used to quantizy the space in the x direction 
%
%  stepY     := the step used to quantizy the space in the y direction
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

if nargin<5 || isempty(xrange)
    xrange=[min(x), max(x)];
end

if nargin<6 || isempty(yrange)
    yrange=[min(y), max(y)];
end

xrange=xrange+[-1 1]; %enlarge the range to keep in the histogram the border like situation
yrange=yrange+[-1 1]; %enlarge the range to keep in the histogram the border like situation

x=x-[repmat(xrange(1),size(x,1),1)];
y=y-[repmat(yrange(1),size(y,1),1)];

m=zeros(nbins,mbins);

stepX=(xrange(2)-xrange(1))/(nbins);
stepY=(yrange(2)-yrange(1))/(mbins);
        
for j=1:numel(x)
    xx=ceil(x(j)/stepX);
    yy=ceil(y(j)/stepY);
    
    if xx==0 
        xx=1;
    end
    
    if xx>nbins
        xx=nbins;
    end
    
    if yy==0 
        yy=1;
    end
    
    if yy>mbins
        yy=mbins;
    end
    
    m(yy,xx)=m(yy,xx)+1;
end

m=m./numel(x);

end

