function value = isdiag (X)
%ISDIAG Check diagonality 
%
%  - Input variable(s) -
%  X: a vector or a matrix
%
%  - Output variable(s) -
%  VALUE: contains a code for the diagonality of X.
%  * VALUE=2 if X is a vector. (it is considered as the diagonal vector of
%    the underlying matrix)
%  * VALUE=1 if X is a diagonal matrix.
%  * VALUE=0 if X is none of the above.
%
%  - Construction -
%  VALUE = ISDIAG(X) returns a code in VALUE for the diagonality of X. 

    N=size(X);
    if (N(1)==1 || N(2)==1) && length(N)<3
        value=2;
	elseif length(N)<3 && isequal(X,diag(diag(X))) 
        value=1;
    else
        value=0;
    end

% faster code when huge > 10000 diagonal matrix
% if nnz(A)>size(A,1)     %fast precheck
% 	result=false;
%     return
% end
% 
% [I,J] = find(A);        
% 
% if ~isempty(I)
%   result = all(I == J);
% else
%   result = true;        %zero matrix is diagonal
% end 