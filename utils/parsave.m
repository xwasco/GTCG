function parsave( varargin)

for i=2:2:numel(varargin)
    a283199231284177asda12389asd23428989asd=varargin{i};
    str = [varargin{i+1},'=a283199231284177asda12389asd23428989asd;'];
    eval(str);
    clear a283199231284177asda12389asd23428989asd;0
    %save('C:\Users\svascon\Dropbox\Matlab\ACCV_Journal\results_Frustum_CVIU\mkmkm.mat',,'-append');
    if exist(varargin{1})
        %str=strrep(['save(''' varargin{1} ''',''' varargin{i+1} ''',''-append'',''-v7.3'');'],'\','\\');
        str=strrep(['save(''' varargin{1} ''',''' varargin{i+1} ''',''-append'');'],'\','\\');
    else
        str=strrep(['save(''' varargin{1} ''',''' varargin{i+1} ''',''-v7.3'');'],'\','\\');    
    end
    
    eval(str);
end

end

