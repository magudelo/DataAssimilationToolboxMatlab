function value = issym(matrix,tol)
%ISSYM Check symmetry 
%
%  - Input variable(s) -
%  MATRIX: Square matrix (mxm) or an array of l square matrices (m x m x l)
%
%  TOL: Tolerance to decide on symmetry. If difference between two matrix 
%  elements is within this tolerance, the elements are considered
%  identical.
%
%  - Output variable(s) - 
%  VALUE: Indicates symmetry, 1=symmetric, 0=not symmetric.
%
%  - Construction -
%  VALUE = ISSYM(MATRIX,TOL) indicates if MATRIX is a symmetric matrix
%  within the tolerance TOL.
%
%  VALUE = ISSYM(MATRIX) indicates if MATRIX is a symmetric matrix.
%  Since the tolerance is not provided, itis determined automatically.
%  


    ni=nargin;

    N=size(matrix);

    if length(N)==2 && N(1)==N(2)       %Square matrix

        if ni<2;tol = 10*eps(max(abs(diag(matrix))));end;           % no tolerance provided: use 10*eps of max diag element
                                                                    % also see comment on bottom        
        value=all(all( abs(matrix - matrix')<=tol ));               % value=1: symmetric within tolerance

    elseif length(N)==3 && N(1)==N(2)   %3D Square matrix array

        if ni==2 && tol==0
            value=isequal(matrix,permute(matrix,[2 1 3]));          % tol=0: must be exactly symmetric
        else
            tols=zeros(1,N(3));
            for i=1:N(3)
                tmpmatrix = matrix(:,:,i);
                if ni<2
                    tol = 10*eps(max(abs(diag(tmpmatrix))));
                end
                tols(i)=all(all( abs(tmpmatrix - tmpmatrix')<=tol ));
            end
            value=all(tols);
        end

    else
        value=0;
    end
    
end

%better would be using eps of max of each off-diagonal entry in a
%for loop but would be to extensive. Diag is made as a choise.
