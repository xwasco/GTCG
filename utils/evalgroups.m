function [Ps,Rs,TPs,FPs,FNs] = evalgroups(group,GTs,crit)
% EVALGROUPS computes precision and recall scores, as well as the number of
% true positives, false positives and false negatives, give a set of
% detected groups and a set of ground truth groups.
%
% Thanks Francesco Setti to provide me this evaluation tool.
% This code has been used in this work to evaluate the results, please cite [18] 
% if you use this code to evaluate your algorithms. 
%
% With respect to the original algorithm now it takes into account multiple
% groundtruth for the same dataset that can be provided by different
% annotators.
%
% [18]Setti, F., Hung, H., Cristani, M.: Group Detection in Still Images by F-formation Modeling: a Comparative Study. In: WIAMIS. (2013) 
% 
% INPUT:
% ======
%  group      := <G-elements> cell array, each element contains a group
%                detected by tested algorithm, each one defined by the
%                array of subjects' ID of individuals belonging to it.
%  GT         := <G-elements> cell array, each element contains a group
%                provided by ground truth, each one defined by the array of
%                subjects' ID of individuals belonging to it. 
%                In case that multiple groundtruth are available (due to
%                multiple annotators) GT is a nxF cell array in which each

%  crit       := string, it defines the criterium to establish whether
%                group has been detected or not. It can be 'ALL' or 'CARD':
%                 - 'CARD' : the detected group have to contain at least
%                            2/3 of the elements of GT group (and
%                            vice-versa).
%                 - 'ALL'  : the detected group and GT group have to match
%                            perfectly.
%                 - 'ZANOTTOBMVC2011': at least the 60% of the member are
%                 correctly groupend (T=0.6);
% 
% OUTPUT:
% =======
% All the output are vectors   
%  Ps         := the ability to detect 'only' the ground truth (or the
%                ability not to generate false positives).
%                   Pr = TP / (TP+FP)
%  Rs         := the ability to detect 'all' the ground truth (or the
%                ability not to generate false negatives).
%                   Re = TP / (TP+FN)
%  TPs        := number of True Positives
%  FPs        := number of False Positives
%  FNs        := number of False Negatives
% 


%% INITIALIZATION
% Default values definition and input data manipulation

% DEFAULT: crit = 'card'
if ~exist('crit','var')
    crit = 'card' ;
elseif ~ischar(crit)
    warning('evalgroups:crit_warning','The input parameter ''crit'' in function evalgroups has to be a string with values ''card'' or ''all''! ''crit'' will be redefined according to the default value')
    crit = 'card' ;
end

% CHECK: group
if isempty(group)
    warning('group_empty','Input variable GROUPS is empty! Precision should be NaN!')
end

% CHECK: GT
if isempty(GTs)
    warning('GT_empty','Input variable GT is empty! Recall should be NaN!')
end


%% PROCESSING 

%
TPs=zeros(size(GTs,1),1);
FPs=zeros(size(GTs,1),1);
FNs=zeros(size(GTs,1),1);
Ps=zeros(size(GTs,1),1);
Rs=zeros(size(GTs,1),1);

for k=1:size(GTs,1)

    GT=GTs{k,:};
    
% Initialize the 'True Positive'
TP = 0 ;

% For each GT group
for ii = 1:size(GT,2)
    
    % Select the GT group and its cardinality (= number of elements)
    GTtmp  = GT{ii} ;
    GTcard = length(GTtmp) ;
    
    % for each detected group
    for jj = 1:numel(group) %size(group,2)
        
        % Select the detected group and its cardinality (= number of elements)
        grouptmp  = group{jj} ;
        groupcard = length(grouptmp) ;
        
        switch crit
            
            case 'half' 
                %condition to align the result to the paper:
                %Discovering Groups of People in Images 
                %Wongun Choi, Yu-Wei Chao, Caroline Pantofaru and Silvio Savarese

                % Find the intersection between the GT and detected groups
                % and its cardinality.
                inters     = intersect(GTtmp,grouptmp) ;
                interscard = length(inters) ;
                % If the cardinality of GT and detected groups are 2, the
                % groups have to match perfectly.
                if groupcard==2 && GTcard == 2
                    if isempty(setxor(GTtmp,grouptmp))
                        TP = TP + 1 ;
                        break
                    end
                elseif ge(interscard/max(GTcard,groupcard),1/2)
                    TP = TP + 1 ;
                    break
                end
                % NOTE: this condition is maybe not needed because it could
                % be implicit in the computation of cardinality condition.
                % To be verified!!!
            
            case 'card'
                % Find the intersection between the GT and detected groups
                % and its cardinality.
                inters     = intersect(GTtmp,grouptmp) ;
                interscard = length(inters) ;
                % If the cardinality of GT and detected groups are 2, the
                % groups have to match perfectly.
                if groupcard==2 && GTcard == 2
                    if isempty(setxor(GTtmp,grouptmp))
                        TP = TP + 1 ;
                        break
                    end
                elseif ge(interscard/max(GTcard,groupcard),2/3-eps)
                    TP = TP + 1 ;
                    break
                end
                % NOTE: this condition is maybe not needed because it could
                % be implicit in the computation of cardinality condition.
                % To be verified!!!
                
            case 'DPMM'
                % Find the intersection between the GT and detected groups
                % and its cardinality.
                inters     = intersect(GTtmp,grouptmp) ;
                interscard = length(inters) ;
                % If the cardinality of GT and detected groups are 2, the
                % groups have to match perfectly.
                if groupcard==2 && GTcard == 2
                    if isempty(setxor(GTtmp,grouptmp))
                        TP = TP + 1 ;
                        break
                    end
                elseif ge(interscard/max(GTcard,groupcard),0.6-eps)
                    TP = TP + 1 ;
                    break
                end
                % NOTE: this condition is maybe not needed because it could
                % be implicit in the computation of cardinality condition.
                % To be verified!!!                
                               
            case 'all'
                if isempty(setxor(GTtmp,grouptmp))
                    TP = TP + 1 ;
                end
                
            otherwise
                error('evalgroups:crit_error','The input parameter ''crit'' in function ff_evalgroups has to be a string with values ''card'' or ''all''!')
        end
    end
end

% Define the number of false positives and negatives.
FP = numel(group)-TP; %size(group,2) - TP ;%
FN = size(GT,2) - TP ;

% Define precision.
precision = TP ./ (TP+FP) ;
if isnan(precision)
    precision=0;
end

% Define recall.
recall = TP ./ (TP+FN) ;
if isnan(recall)
    recall=0;
end
% NOTE: this control on isnan is needed in order to avoid problems in
% f-formation main code, but it should be changed in order to take into
% account that this 0 is different from the one generated by TP=0 and
% FP~=0 or FN~=0.

TPs(k)=TP;
FPs(k)=FP;
FNs(k)=FN;
Ps(k)=precision;
Rs(k)=recall;

end
