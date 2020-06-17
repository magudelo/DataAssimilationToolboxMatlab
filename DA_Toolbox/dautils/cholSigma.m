function [U,code] = cholSigma(Sigma,reqSemiDef,isDiag,lower)
%CHOLSIGMA Calculate Cholesky decomposition 
%
%  - Input variable(s) -
%  SIGMA: can be either a square symmetric matrix or a vector. When SIGMA 
%  is a vector, it contains the diagonal elements of a diagonal matrix.
%
%  REQSEMIDEF: when equal to 1 indicates if positive semi-definiteness of
%  SIGMA should also be checked.
%
%  ISDIAG: To increase performance, it is possible to indicate diagonality
%  of Sigma with isDiag. If isDiag=1 SIGMA is a diagonal matrix, if isDiag=2
%  SIGMA is a diagonal vector, if isDiag=0 SIGMA is a regular matrix.
%  
%  LOWER: If 1 returns the lower triangular matrix; if 0 the upper
%  triangular matrix. (Default =0)
%
%  - Output variable(s) -
%  U: the upper/lower triangular Cholesky decomposition of SIGMA. Is empty 
%  ([ ]) if SIGMA is not positive definite and REQSEMIDEF is 0, or if SIGMA 
%  is indefinite and REQSEMIDEF is 1.
%
%  CODE: indicates the kind of definiteness of SIGMA. Namely:
%  * 0: Positive definite
%  * 2: Not positive definite (U=[])
%  * 1: Positive semi-definite (reqSemiDef=1)
%  * 3: Indefinite (reqSemiDef=1, U=[])
%
%  - Construction -
%  [U,code] = CHOLSIGMA(SIGMA, reqSemiDef, isDiag, lower) returns the upper 
%  or lower triangular Cholesky decomposition U of SIGMA.
%
%  [U,code] = CHOLSIGMA(SIGMA, reqSemiDef, isDiag) returns the upper 
%  triangular Cholesky decomposition U of SIGMA.
%
%  [U,code] = CHOLSIGMA(SIGMA,reqSemiDef) returns the upper triangular 
%  Cholesky decomposition U without using the knowledge of isDiag.
%
%  [U,code] = CHOLSIGMA(SIGMA) returns the upper triangular Cholesky
%  decomposition without using the knowledge of isDiag and will not calculate
%  the decomposition in case SIGMA is a positive semi-definite matrix.
%

    ni=nargin;
    if ni<2;reqSemiDef=1;end;
    if ni<3;isDiag=0;end;
    if ni<4;lower=0;end;

    retMat=0;
    %Diagonal matrix supplied: change into vector
    if isDiag==1 || (min(size(Sigma))~=1 && isequal(Sigma,diag(diag(Sigma))) )
       Sigma=diag(Sigma); 
       isDiag=2;
       retMat=1;
    end
    
    %Diagonal matrix or vector supplied
    if isDiag==2 || size(Sigma,1)==1 || size(Sigma,2)==1   
        
        tol=10*eps(max(abs(Sigma)));
        negs = sum(Sigma < -tol);
        zeros = sum(Sigma < tol);
        
        if (negs==0) && (zeros==0)
            code = 0;   % Positive definite
            U=sqrt(Sigma);
            if retMat;U=diag(U);end;
            
        elseif reqSemiDef
            
            if (negs==0) && (zeros~=0)
                code = 1;   % Positive semi-definite
                U=sqrt(Sigma);
                if retMat;U=diag(U);end;
            else
                code = 3;   % Indefinite
                U = [];
            end
            
        else
            code = 2; % Not positive definite
            U = [];
        end
    %Non-Diagonal matrix: matrix supplied
    else        
    
        if size(Sigma,1)~= size(Sigma,2)
            code=10; U=[]; % Not square
            return;
        end

        if ~issym(Sigma)
            code=11; U=[]; % Not symmetric
            return;
        end

        if lower 
            [U,err] = chol(Sigma,'lower');
        else
            [U,err] = chol(Sigma);
        end

        if err > 0  % Not positive definite: check semi definite
                   
            if reqSemiDef
                [V,D] = eig(Sigma);     %Sigma= V D V'

                D = diag(D);
                tol = eps(max(D)) * length(D);
                ind = (abs(D) > tol);
                D(~ind) = 0;             % set small eigenvalues=0
                negs = sum(D<0);        % number of negative eigenvalues 

                if (negs==0)            %no negative eigenvalues
                    if lower
                        U = V*diag(sqrt(D));
                    else
                        U = diag(sqrt(D))*V';
                    end
                    code = 1;   % Positive semi-definite
                else
                    code = 3;   % Indefinite
                    U = [];
                end
            else
                U = [];
                code = 2; % Not positive definite
            end
        else
            code = 0; % Positive definite
        end
    end
end

