function [matrix] = interpolate3D(array3D,index,reqVal,varargin)
%INTERPOLATE3D interpolates to find MATRIX, the values of ARRAY3D 
% at the point REQVAL in the INDEX vector.
%
%  - Input variable(s) -
%  ARRAY3D: a 3D array
%  INDEX: the interpolation index for REQVAL
%  REQVAL: the desired value in INDEX
%  EXTRAPOL: 0=interpolate, 1=extrapolate
%
%  - Output variable(s) -
%  MATRIX: the interpolated matrix
%
%  - Construction -
%  [MATRIX] = INTERPOLATE3D(ARRAY3D,INDEX,REQVAL,EXTRAPOL) interpolates or 
%  extrapolates to find MATRIX, the values of ARRAY3D at the point REQVAL
%  in the INDEX vector.
%
%  [MATRIX] = INTERPOLATE3D(ARRAY3D,INDEX,REQVAL) interpolates to 
%  find MATRIX, the values of ARRAY3D at the point REQVAL in the INDEX vector.

% Note:Matlab uses copy-on-write to avoid making a copy of the input argument inside 
% the function workspace until or unless you modify the input argument! So
% no problem to use complete 3D array.

narginchk(3, 4);
ni=nargin;

if isempty(array3D)|| isempty(index)|| isempty(reqVal)
    error('DA:dautils:findHighest:interpolate3D:parMismatch','Parameter missing.');
end

if ni==3;extrapolate=0;end;
if ni==4
    if varargin{1}==1
        extrapolate=1;
    elseif varargin{1}==1
        extrapolate=0;
    else
        error('DA:dautils:findHighest:interpolate3D:parMismatch','Wrong parameter value for interpolate argument.');
    end
end

[valLow,indLow] = findLowest(index,reqVal);
[valHigh,indHigh] = findHighest(index,reqVal);

if reqVal<index(1)                  %desiredVal too low
    if extrapolate==1
        matrixLow=array3D(:,:,1);
        matrixHigh=array3D(:,:,2);
        valLow=index(1);
        valHigh=index(2);
        matrix = matrixLow+(matrixHigh-matrixLow).*((reqVal-valLow)./(valHigh-valLow));
    else
        matrix=array3D(:,:,1);
    end
elseif reqVal>index(end)            %desiredVal too high
    if extrapolate==1
        matrixLow=array3D(:,:,end-1);
        matrixHigh=array3D(:,:,end);
        valLow=index(end-1);
        valHigh=index(end);
        matrix = matrixHigh+(matrixHigh-matrixLow).*((reqVal-valHigh)./(valHigh-valLow));
    else
        matrix=array3D(:,:,end);
    end    
elseif indLow==indHigh              %right on an index value
    matrix=array3D(:,:,indHigh);
else                                %in between: interpolate
	matrixLow=array3D(:,:,indLow);
    matrixHigh=array3D(:,:,indHigh);
    matrix = matrixLow+(matrixHigh-matrixLow).*((reqVal-valLow)./(valHigh-valLow));
end


