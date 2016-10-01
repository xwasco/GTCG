function dist=JSDiv(P,Q)
% Jensen-Shannon divergence of two probability distributions
%  dist = JSD(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each
% P =  1 x nbins
% Q =  1 x nbins
% dist = n x 1


if size(P,2)~=size(Q,2)
    error('JSDivergence: the number of columns in P and Q should be the same\n');
end

% normalizing the P and Q
Q = Q ./sum(Q);
P = P ./sum(P);

%T=(P>0)&(Q>0);
%if sum(T(:))>0
    M = 0.5.*(P + Q);
    
    dist = 0.5.*KLDiv(P,M) + 0.5*KLDiv(Q,M);
%else
%    dist=inf;
%end
end