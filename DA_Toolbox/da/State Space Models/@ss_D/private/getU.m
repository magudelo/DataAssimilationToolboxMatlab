function newU=getU(u,index)
%GETU obtains input column vector from matrix
%
%  - Input variable(s) -
%  U: matrix of input column vectors
%
%  INDEX: selected column number
%
%  - Output variable(s) -
%  NEWU: selected input column vector
%  
%  - Construction -          
%  NEWU=GETU(U,INDEX) obtains the input column vector defined by INDEX value 
%  from matrix U. If U is empty returns empty matrix [].

	if isempty(u)
        newU=[];
    else
        newU=u(:,index);
	end
            
end