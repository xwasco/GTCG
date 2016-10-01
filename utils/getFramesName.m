function fNames=getFramesName(f)
% GETFRAMESNAME Return a list of file names for a given dataset 
%corresponding to the frames

    listing = dir(f);

    fNames=[];
    for i=1:numel(listing)
        if strcmp(listing(i).name,'.')~=1 && strcmp(listing(i).name,'..')~=1
            fNames=[fNames ;  [i,{listing(i).name}]];
        end
    end

end