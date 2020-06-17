function [nearVal, index] = findNearest(x, desiredVal)
%FINDNEAREST find the closest value to DESIREDVAL that exists in X
%
%  - Input variable(s) -
%  X: a vector or matrix
%  DESIREDVAL: the desired value in x
%
%  - Output variable(s) -
%  LOWVAL: the nearest value
%  INDEX: the index where the nearest value is found
%
%  - Construction -
%  [LOWVAL, INDEX] = FINDHIGHEST(X, DESIREDVAL) finds the nearest value
%  to DESIREDVAL that exsists in X and returns both value and index. 
 
x = x(:)'; %% this resizes the matrix
if nargin == 2
    if ismember(x,desiredVal) == 1 %% This is the O(1) case
        index = find(x == desiredVal);
        nearVal = x(x==desiredVal);
    else 
        [~, index] = min(abs(desiredVal-x)); 
        nearVal = x(index);
    end
else
    error('DA:dautils:findNearest:argMismatch','You have not entered the correct amount of parameters');
end