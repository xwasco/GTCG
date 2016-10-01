function [ H ] = weights_entropy( x )
%ENTROPY Given a probability distribution x compute the entropy

lx=log2(x+eps);

H=-sum(x.*lx);

end

