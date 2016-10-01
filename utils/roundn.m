function [ val ] = roundn( val, ndec )
%ROUNDN Round a number to n-decimals.

if nargin<2
    ndec = 2;
end

f = 10.^ndec;
val = round(f*val)/f; 

end

