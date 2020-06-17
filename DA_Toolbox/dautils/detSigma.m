function value = detSigma(Sigma)
%DETSIGMA Calculate determinant of covariance matrix Sigma.
%
%  - Input variable(s) -
%  SIGMA: a matrix or a vector that contains the diagonal elements of a 
%  diagonal matrix.
%
%  - Output variable(s) -
%  VALUE: the value of the determinant of SIGMA.
%
%  - Construction -
%  DETSIGMA(SIGMA) returns the determinant of SIGMA in a more efficient 
%  manner than the function det(sigma).

    if min(size(Sigma))==1                  %vector provided: assume diagonal matrix
        value=prod(Sigma);
    elseif isequal(Sigma,diag(diag(Sigma)))	%diagonal matrix
        value=prod(diag(Sigma));
    elseif isequal(triu(Sigma),Sigma)    	%upper triangular matrix
        value=prod(diag(Sigma));
    elseif isequal(tril(Sigma),Sigma)       %lower triangular matrix
        value=prod(diag(Sigma));        
    else
        value=det(Sigma);
    end
    
end