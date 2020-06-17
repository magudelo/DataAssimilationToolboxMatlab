function array=findArrayVal(array3D,tIndex,tMethod,t)
%FINDARRAYVAL Find the correct 2D array from a 3D array.
%
%  - Input variable(s) -
%  ARRAY3D: a 3D array where the third dimension represents different times.
%
%  T: requested time value.
%
%  TINDEX: time index vector used to indicate which time corresponds to
%  each vector or matrix in the 3D arrays.
%
%  TMETHOD: string that indicates the rounding/interpolation
%  method to retrieve the covariance matrix or mean vector from the
%  3D array at time T while taking into account tIndex. The
%  following strings are possible:
%  1) 'low': find 2D array corresponding to the lower bound time in tIndex 
%  2) 'high': find 2D array corresponding to the higher bound time in tIndex
%  3) 'near': find 2D array corresponding to the nearest time in tIndex
%  4) 'interp': find interpolated 2D array corresponding to tIndex. When
%     T is lower/higher than lowest/highest value in tIndex, the first/last
%     element is used.
%  5) 'extrap': find interpolated 2D array corresponding to tIndex. When
%     T is lower/higher than lowest/highest value in tIndex, extrapolation
%     is applied.
%
%  - Output variable(s) -
%  ARRAY: the 2D array selected from the 3D array that corresponds to time T.
%
%  - Construction -
%  ARRAY=FINDARRAYVAL(ARRAY3D,TINDEX,TMETHOD,T) returns a 2D array from the
%  3D array corresponding to time T. To locate the correct 2d array, its index
%  is set identical to the index of T in TINDEX. When T is not exactly found
%  in TINDEX, the method defined in TMETHOD is used to find the correct
%  index.
%
%  For example:
%     If ARRAY3D is the 2 x 2 x 3 array of [1 2;2 1], [2 3;3 2] and 
%     [3 4;4 3] and TINDEX=[2 3 5] and T=2.4 than TMETHOD
%     'low'     yields [1 2;2 1]
%     'high'    yields [2 3;3 2]
%     'near'    yields [1 2;2 1]
%     'interp'  yields [1.4 2.4;2.4 1.4]
%     'extrap'  yields [1.4 2.4;2.4 1.4]
%     when T=1, TMETHOD
%     'interp'  yields [1 2;2 1]
%     'extrap'  yields [0.6 1.6;1.6 0.6]

	narginchk(4, 4);

    switch tMethod
        case 'low'
            [~,index]=findLowest(tIndex, t);   
            array=array3D(:,:,index);
        case 'high'
            [~,index]=findHighest(tIndex, t);     
            array=array3D(:,:,index);            
        case 'near'
        	[~,index]=findNearest(tIndex, t);  
            array=array3D(:,:,index);            
        case 'interp'
            array=interpolate3D(array3D,tIndex,t);
        case 'extrap'
            array=interpolate3D(array3D,tIndex,t,1);
    end
    
end        