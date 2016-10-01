function dist=KLDiv(P,Q)
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1

Q = Q ./sum(Q);

P = P ./sum(P); %repmat(sum(P,2),[1 size(P,2)]);

dt=log2(Q);
d1=(P.*dt);
d1(isnan(d1))=0;
d1=-sum(d1);

d2=P.*log2(P);
d2(isnan(d2))=0;
d2=sum(d2);

dist=d1+d2;

end

