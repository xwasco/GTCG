function plotresult(PP,PPh,DD,c,ii,nn)
% PLOTRESULT plots part of the results of BENSOLVE, the upperimage of (P)
% and the lower image of (D^*) for 2 or 3 objectives
%
% USAGE:
%
% plotresult(PP,PPh,DD)
% plotresult(PP,PPh,DD,c,ii,nn)
%
% The geometric duality parameter 'c' is required just in case it was
% specified in BENSOLVE, i.e. if 'c' apart from its default value
% c=(1,...,1)^T was used.
%
% 'ii' and 'nn' are optional arguments to specify subplots, see example4.m
%
% LICENSE: This file is part of BENSOLVE, see bensolve.m

    display('Plotting the result of BENSOLVE ...')

    eps=1e-10;    
    q=size(PP,1);
    
    if ~exist('c','var') || isempty(c)
        c=ones(q,1);
    end
    
    if ~exist('ii','var')
        ii=1;
    end    
    if ~exist('nn','var')
        nn=1;
    end
    
    %%% plot upper image of (P)
    subplot(nn,2,2*ii-1);
   
    d_hat=sum(DD(1:end-1,:),2)./size(DD,2);            
    eta=[d_hat;1 - d_hat' * c(1:q-1,1) ];    
    
    if q > 3 
        error('Dimension too large for plot')
    elseif q < 2 
        error('No plot of one-dimensional or empty objects')
    elseif q==3  
        d=0;
        for i=1:size(PP,2)
            for j=1:size(PP,2)
                dist=norm(PP(:,i)-PP(:,j));
                if d < dist
                    d = dist;
                end
            end
        end 
         if d>eps
            gamma=max(eta'*PP)+(2/3*d)*min(eta'*PPh);
        else
            gamma=max(eta'*PP)+1;
        end                
        ZZ=PP;
        for i=1:size(PP,2)
            for j=1:size(PPh,2)
                ZZ=[ZZ,PP(:,i) + (gamma - eta'*PP(:,i))/(eta'*PPh(:,j))*PPh(:,j)]; %#ok
            end
        end  
        Kp=convhulln(ZZ'); 
        x=ZZ(1,:)';
        y=ZZ(2,:)';
        z=ZZ(3,:)';
        warning off %#ok
        tr=TriRep(Kp,x,y,z);
        warning on %#ok
        fe = featureEdges(tr,pi/50)';
        trisurf(tr, 'FaceColor', 'cyan', 'EdgeColor','none', 'FaceAlpha', 0.8); axis equal;
        hold on; 
        plot3(x(fe), y(fe), z(fe), 'k', 'LineWidth',1.5); 
        scatter3(PP(1,:),PP(2,:),PP(3,:),'MarkerEdgeColor','black','MarkerfaceColor','black');                
        hold off;
        view([-20 -10])
    elseif q==2
        d=0;
        for i=1:size(PP,2)
            for j=1:size(PP,2)
                dist=norm(PP(:,i)-PP(:,j));
                if d < dist
                    d = dist;
                end
            end
        end
        if d>eps
            gamma=max(eta'*PP)+(2/3*d)*min(eta'*PPh);
        else
            gamma=max(eta'*PP)+1;
        end
        ZZ=PP;
        for i=1:size(PP,2)
            for j=1:size(PPh,2)
                ZZ=[ZZ,PP(:,i) + (gamma - eta'*PP(:,i))/(eta'*PPh(:,j))*PPh(:,j)]; %#ok
            end
        end   
        x=ZZ(1,:)';
        y=ZZ(2,:)';  
        id=convhull(x,y);
        fill(x(id),y(id),'cyan');
        axis auto;
        hold on;
        scatter(PP(1,:),PP(2,:),'MarkerEdgeColor','black','MarkerfaceColor','black');
        hold off;
    end
    
    title('upper image of (P)');
    
    
    %%% plot lower image of (D^*)
    subplot(nn,2,2*ii);
    
    q=size(DD,1);
    if q > 3 
        error('Dimension too large for plot')
    elseif q < 2 
        error('No plot of one-dimensional or empty set')
    elseif q==3  
        d=max(DD(3,:))-min(DD(3,:));
        if d>eps
            gamma=min(DD(3,:))-(1/2*d);
        else
            gamma=min(DD(3,:))-1;
        end                
        ZZ=DD;
        for i=1:size(DD,2)
            ZZ=[ZZ,[DD(1:2,i);gamma]]; %#ok
        end                 
        Kp=convhulln(ZZ'); 
        x=ZZ(1,:)';
        y=ZZ(2,:)';
        z=ZZ(3,:)';
        warning off %#ok
        tr=TriRep(Kp,x,y,z);
        warning on %#ok
        fe = featureEdges(tr,pi/50)';
        trisurf(tr, 'FaceColor', 'yellow', 'EdgeColor','none', 'FaceAlpha', 0.8); axis equal;
        hold on; 
        plot3(x(fe), y(fe), z(fe), 'k', 'LineWidth',1.5); 
        scatter3(DD(1,:),DD(2,:),DD(3,:),'MarkerEdgeColor','black','MarkerfaceColor','black');                
        hold off;
        view([15 80])
    elseif q==2
        d=max(DD(2,:))-min(DD(2,:));
        if d>eps
            gamma=min(DD(2,:))-(1/4*d);
        else
            gamma=min(DD(2,:))-1;
        end
        ZZ=DD;
        for i=1:size(DD,2)
            ZZ=[ZZ,[DD(1,i);gamma]]; %#ok
        end
        x=ZZ(1,:)';
        y=ZZ(2,:)';  
        id=convhull(x,y);
        fill(x(id),y(id),'yellow');
        axis auto;
        hold on;
        scatter(DD(1,:),DD(2,:),'MarkerEdgeColor','black','MarkerfaceColor','black');        
        hold off;        
    end
    
    title('lower image of (D^*)');    
    
    display('done.')
end


