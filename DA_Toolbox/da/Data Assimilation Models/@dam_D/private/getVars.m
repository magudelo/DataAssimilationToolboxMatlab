function vars=getVars(cov,full,ind,dataind)
%GETVARS obtains variances from variance matrix or covariance 3D array
%
%  - Input variable(s) -
%  COV: variance matrix or covariance 3D array
%  FULL: inidication variance matrix (=0) or covariance 3D array (=1)
%  IND: index of required states
%  DATAINDEX: index of required points
%
%  - Output variable(s) -
%  VARS: matrix of variances
%  
%  - Construction -          
%  VARS=GETVARS(COV,FULL,IND) obtains variances from variance matrix 
%                         or covariance 3D array as specified in IND

        if full
            vars=diag3D(cov);
            vars=vars(ind,dataind);
        else
            vars=cov(ind,dataind);
        end
            
end