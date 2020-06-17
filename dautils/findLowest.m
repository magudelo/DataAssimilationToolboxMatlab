function [lowVal, index] = findLowest(x, desiredVal)
%FINDLOWEST find the lower bound value to DESIREDVAL that exists in X
%
%  - Input variable(s) -
%  X: a vector or matrix
%  DESIREDVAL: the desired value in x
%
%  - Output variable(s) -
%  LOWVAL: the lower bound value
%  INDEX: the index where the lower bound value is found
%
%  - Construction -
%  [LOWVAL, INDEX] = FINDLOWEST(X, DESIREDVAL) finds the lower bound value 
%  to DESIREDVAL that exsists in X and returns both value and index. 

x = x(:)'; %% this resizes the matrix
if nargin == 2
    index = find(x == desiredVal,1); 
    if ~isempty(index)
        lowVal = x(x==desiredVal);
    else 
        index=find(x>=desiredVal,1)-1;
        if isempty(index)
            index = length(x);
        elseif index<1
            index=1;
        end
        lowVal=x(index);
    end
else
    error('DA:dautils:findLowest:argMismatch','You have not entered the correct amount of parameters');
end