function diags = diag3D(x)
%DIAG3D Returns the diagonal elements of 3D array of matrices.
%
%  - Input variable(s) -
%  X: a 3D array of matrices
%
%  - Output variable(s) -
%  DIAGS: the diagonal elements. Each column represents the diagonals of 
%  one matrix from the 3D array.
%
%  - Construction -
%  DIAGS = DIAG3D(X) returns the diagonal elements of 3D array of matrices.

    z = size(x);
    m = z(1);
    n = z(2);
    mn = m * n;
    p = prod(z)/mn;
    k = (1:m+1:m*n)';
    k = repmat(k,1,p);
    k = bsxfun(@plus,k,(0:(p-1))*mn);
    diags = x(k);
end