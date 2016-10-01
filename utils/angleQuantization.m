function [ v ] = angleQuantization( theta, range,numBin )
    R=linspace(range(1),range(2),numBin);
    h=hist(theta,R);
    v=R(find(h==1));
end

