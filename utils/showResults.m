function showResults( precisions,recalls,id )
%SHOWRESULTS Show the results

if nargin<3
    id=0;
end

for i=1:size(precisions,2)
    fprintf(['Groundtruth #' num2str(i+id) ': ']);
    fprintf('Pre: %d ',mean(precisions(:,i)));
    fprintf('Rec: %d ',mean(recalls(:,i)));
    f1=2*(mean(precisions(:,i)).*mean(recalls(:,i)))./(mean(precisions(:,i))+mean(recalls(:,i)));    
    fprintf('F1: %d\n',f1);
end

if size(precisions,2)>1
    mprec=mean(reshape(precisions,1,numel(precisions)));
    mrec=mean(reshape(recalls,1,numel(recalls)));
    f1=2*(mprec*mrec)./(mprec+mrec);   
    fprintf(['Average over ' num2str(size(precisions,2)) ' groundtruths Pre: ' num2str(mprec) ' Rec: ' num2str(mrec) ' F1: ' num2str(f1) '\n']);
end

end

