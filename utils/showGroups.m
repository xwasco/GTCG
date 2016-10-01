function showGroups(fid,frames,detection,gt,param)

    figure(param.show.groups);
    
    if isfield(param, 'homography')
        for i=1:size(frames,1)
            t=[fs(i,2:3) , 1];
            t=t*T;
            FL=[fs(i,2)+param.frustum.length*cos(fs(i,4)),fs(i,3)+param.frustum.length*sin(fs(i,4)) , 1]*T;
            FL=FL/FL(3);
            fs(i,2:3)=t(1:2)/t(3);
            fs(i,5:6)=FL(1:2);
        end
    else
        fs(:,5:6)=[fs(:,2)+param.frustum.length.*cos(fs(:,4)),fs(:,3)+param.frustum.length.*sin(fs(:,4))];
    end

    for i=1:numel(fid)
        f=fid(i);
        if isfield(param, 'framesNames') && ~isempty(param.framesNames)
            %display the image of the i-th frame
            imsow(imread([param.datasetDir '/' param.framesDir '/' param.framesNames{f,2}]));
            hold on
        end
    end
end