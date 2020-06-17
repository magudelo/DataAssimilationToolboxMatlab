function [highVal, index] = findHighest(x, desiredVal)
%FINDHIGHEST find the higher bound value to DESIREDVAL that exists in X
%
%  - Input variable(s) -
%  X: a vector or matrix
%  DESIREDVAL: the desired value in x
%
%  - Output variable(s) -
%  LOWVAL: the higher bound value
%  INDEX: the index where the higher bound value is found
%
%  - Construction -
%  [LOWVAL, INDEX] = FINDHIGHEST(X, DESIREDVAL) finds the higher bound value 
%  to DESIREDVAL that exsists in X and returns both value and index. 

x = x(:)'; %% this resizes the matrix
if nargin == 2
    index = find(x == desiredVal,1); 
    if ~isempty(index)
        highVal = x(x==desiredVal);
    else 
        index=find(x>=desiredVal,1);
        if isempty(index)
            index = length(x);
        elseif index<1
            index=1;
        end
        highVal=x(index);
    end
else
    error('DA:dautils:findHighest:argMismatch','You have not entered the correct amount of parameters');
end