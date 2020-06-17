classdef nm_gauss_ltv < nm_gauss
    %NM_GAUSS_LTV Constructs a multivariate Linear Time Variant Gaussian 
    %Noise Model object
    %
    %  NM_GAUSS_LTV Properties:
    %
    %     p - amount of variables (read only)  
    %
    %     l - length of 3D-array(s), the third dimension (read only)   
    %
    %     mu - 3D array of mean column vectors where the third dimension
    %     represents step numbers. For example: for a bivariate distribution 
    %     with array length 3 mu is a 2x1x3 array (PxLxL). If mu is LTI, a 
    %     column vector as in NM_GAUSS_LTI can be used. 
    %
    %     Sigma - 3D array of covariance matrices where the third dimension
    %     represents step numbers. (Symmetric matrices or column vectors 
    %     that represent the diagonal elements of diagonal matrices.)
    %     For example: for a bivariate distribution with array 
    %     length 3, Sigma is a 2x2x3 array (PxPxL), or a 2x1x3 array (PxLxL)
    %     in case of diagonal vectors. 
    %     If Sigma is LTI, a matrix or column vector as in NM_GAUSS_LTI can
    %     be used. 
    %
    %     kIndex - step number vector used to indicate which step number  
    %     corresponds to each vector or matrix in the 3D arrays. 
    %     (Default = values '0' to 'L-1').
    %     For example:
    %     If SIGMA is the 2 x 2 x 3 array of [1 2;2 1], [2 3;3 2] and 
    %     [3 4;4 3] and kIndex=[2 3 5] then the matrix [2 3;3 2] is the
    %     the covariance matrix at k=3.
    %     
    %     kMethod - string that indicates the rounding
    %     method to retrieve the covariance matrix or mean vector from the
    %     3D array at a step k while taking into account KINDEX. The
    %     following strings are possible (default='low'):
    %     1) 'low': find 2D array corresponding to the lower bound step 
    %               number in KINDEX 
    %     2) 'high': find 2D array corresponding to the higher bound step 
    %               number in KINDEX 
    %     3) 'near': find 2D array corresponding to the nearest step 
    %               number in KINDEX 
    %       For example:
    %       Given the 2 x 2 x 3 array of [1 2;2 1], [2 3;3 2] and 
    %       [3 4;4 3] and KINDEX=[2 3 7] and k=4 than KMETHOD
    %       'low'     yields [2 3;3 2]
    %       'high'    yields [3 4;4 3]
    %       'near'    yields [2 3;3 2] 
    %
    %  NM_GAUSS_LTV Construction:
    %     OBJ = NM_GAUSS_LTV(MU,SIGMA,KINDEX,KMETHOD) creates an object OBJ 
    %     representing a multivariate linear time variant gaussian noise   
    %     model with mean MU and covariance matrix SIGMA. MU or SIGMA must
    %     be a LTV 3D array. If not, use NM_GAUSS_LTI which is more efficient.
    %     If MU and SIGMA are both LTV, the length of their 3D arrays, L, 
    %     must be equal.
    %     
    %     OBJ = NM_GAUSS_LTV(MU,SIGMA,KINDEX) uses 'low' for KMETHOD.
    %
    %     OBJ = NM_GAUSS_LTV(MU,SIGMA,KMETHOD) uses values '0' to 'L-1' for 
    %     KINDEX.    
    %
    %     OBJ = NM_GAUSS_LTV(MU,SIGMA) uses values '0' to 'L-1' for 
    %     KINDEX and 'low' for KMETHOD.      
    %
    %     OBJ = NM_GAUSS_LTV(SIGMA,KMETHOD) creates a multivariate LTV 
    %     gaussian noise model with zero mean and covariance matrix SIGMA.
    %     SIGMA must be a 3D array of symmetric matrices in this case.  
    %     Values '0' to 'L-1' are used for KINDEX.
    %
    %     OBJ = NM_GAUSS_LTV(MU,KMETHOD) creates a multivariate LTV 
    %     gaussian noise model with mean MU and the identity matrix as
    %     covariance matrix. MU must be a 3D array in this case.    
    %     Values '0' to 'L-1' are used for KINDEX.    
    %
    %     OBJ = NM_GAUSS_LTV(SIGMA) creates a multivariate LTV gaussian
    %      noise model with zero mean and covariance matrix SIGMA.
    %     SIGMA must be a 3D array of symmetric matrices in this case.  
    %     Values '0' to 'L-1' are used for KINDEX and 'low' for KMETHOD.
    %
    %     OBJ = NM_GAUSS_LTV(MU) creates a multivariate LTV gaussian
    %     noise model with mean MU and the identity matrix as covariance
    %     matrix. MU must be a 3D array in this case.    
    %     Values '0' to 'L-1' are used for KINDEX and 'low' for KMETHOD.    
    %     
    %     OBJ = NM_GAUSS_LTV(CELL), with CELL a cell array, creates a 
    %     LTV gaussian noise model based on the information in CELL.
    %     
    %  NM_GAUSS_LTV Methods:
    %     cov - Returns the covariance matrix.
    %     var - Returns the variance matrix.
    %     mean - Returns the mean vector.
    %     sample - Draw sample(s) from distribution.
    %     pdf - Returns density from distribution.
    %     cholcov - Returns the upper triangular square root matrix of a 
    %               Positive-Semidefinite Sigma.      
    %
    %  See also NM_GAUSS_LTI, NM_GAUSS_HANDLE    
    
    %mu could have been choosen as matrix, but is also a 3d array for 
    %consistency and to make the user input more error proove: 
    %easier to check if sizes are correct of supplied mu
    %matrix, transpose it when necessary,...
    %Note: memory usage is the same for matrix or 3d vector array
    
    properties (Access = public)
        mu;             %Mean vector
        Sigma;          %Covariance Matrix
        kIndex;         %Step index array
        kMethod='low';  %Step index interpolation/rounding method
    end
    
    properties (Dependent = true, SetAccess = private)
        p;              %Amount of variables 
        l;              %length of 3D-array
    end
    
    methods 
    
        function obj = nm_gauss_ltv(varargin)
        % NM_GAUSS_LTV Constructor    
            narginchk(1, 4);
            ni = nargin;

            if ni==1 && isa(varargin{1},'cell')
                varargin=varargin{1};
                if length(size(varargin))> 2 || (size(varargin,1)>1 && size(varargin,2)>1) 
                    error('DA:NoiseModels:nm_gauss_ltv:cellDim','Cell must be vector sized.')
                else
                    ni=max(size(varargin,1),size(varargin,2));
                end
            end            
            
            dataCnt=0;          %Count the amount of input args of type double 
            for ct=1:ni
                if isa(varargin{ct},'double')
                    dataCnt=dataCnt+1;
                else
                    strCnt=1;
                end
            end
            
            switch ni
                case 1                
                    if dataCnt==1
                        
                        if size(varargin{1},1) == size(varargin{1},2)   % if square matrix: Sigma
                            obj.Sigma=varargin{1};
                            obj.mu=zeros(obj.p,1);
                        else                                            % else: mu
                            obj.mu=varargin{1};
                            obj.Sigma=eye(obj.p);
                        end                        
                        obj.kIndex=(0:obj.l-1);
                        obj.kMethod='low';
                        
                    else
                        error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:argMismatch','Wrong input argument.');
                    end
                case 2
                    if dataCnt==2
                        obj.Sigma=varargin{2};
                        obj.mu=varargin{1};
                        obj.kIndex=(0:obj.l-1);
                        obj.kMethod='low';
                    elseif dataCnt==1 && strCnt==1
                        if size(varargin{1},1) == size(varargin{1},2)   % if square matrix: Sigma
                            obj.Sigma=varargin{1};
                            obj.mu=zeros(obj.p,1);
                        else                                            % else: mu
                            obj.mu=varargin{1};
                            obj.Sigma=eye(obj.p);
                        end     
                        obj.kIndex=(0:obj.l-1);
                        obj.kMethod=varargin{2};                       
                    else
                        error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:argMismatch','Wrong input argument.');
                    end
                case 3
                    if dataCnt==3
                        obj.Sigma=varargin{2};
                        obj.mu=varargin{1};
                        obj.kIndex=varargin{3};
                        obj.kMethod='low';
                    elseif dataCnt==2 && strCnt==1
                        obj.Sigma=varargin{2};
                        obj.mu=varargin{1};
                        obj.kIndex=(0:obj.l-1);
                        obj.kMethod=varargin{3};                       
                    else
                        error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:argMismatch','Wrong input argument.');
                    end                    
                case 4
                    obj.Sigma=varargin{2};
                    obj.mu=varargin{1};
                    obj.kIndex=varargin{3};
                    obj.kMethod=varargin{4};
            end
            
        end        
        
        function obj = set.mu(obj,value)
            
            value=checkArgs(value,'mu');            %check mu                    
            obj.mu=value;
            obj=checkConsist(obj);
            
        end              
        
        function obj = set.Sigma(obj,value)
                                    
            value=checkArgs(value,'Sigma');         %check Sigma                    
            obj.Sigma=value;
            obj=checkConsist(obj);
            
        end            

        function obj = set.kIndex(obj,value)
                                    
            value=checkArgs(value,'kIndex');       %check kIndex                    
            obj.kIndex=value;
            obj=checkConsist(obj);
            
        end            

        function obj = set.kMethod(obj,value)
                                    
            value=checkArgs(value,'kMethod');       %check kMethod                    
            obj.kMethod=value;
            obj=checkConsist(obj);
            
        end                 
        
        function value = get.p(obj)                 %amount of variables 
            
            if isempty(obj.Sigma)
                value = max([size(obj.mu,1) size(obj.mu,2)]);
            else
                value = max([size(obj.Sigma,1) size(obj.Sigma,2)]);
            end
            
        end
        
        function value = get.l(obj)                 %length of 3D array 
            
            value = max([size(obj.Sigma,3) size(obj.mu,3)]);
            
        end
          
        function display(obj)
        %DISPLAY Displays the object summary.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_ltv object     
        %
        %  - Construction -
        %  DISPLAY(OBJ) displays the nm_gauss_ltv object summary.
        
        	disp(' ');
            disp([inputname(1),' = '])
            disp(' ');            
            if ~isempty(obj)             
                disp('Gaussian Noise Model - Time Variant.')
                disp('------------------------------------')
                disp('Properties:');
                if (obj.p<5 && obj.l<3)
                    fprintf('Mean \t\t(.mu) \t = \n');
                    disp(obj.mu);
                    fprintf('Covariance \t(.Sigma) = \n');
                    disp(obj.Sigma);
                    fprintf('Step Indexes (.kIndex) = \n');
                    disp(obj.kIndex);
                else
                    fprintf('Mean \t\t(.mu) \t = \t[%i x %i x %i %s] \n',size(obj.mu,1),size(obj.mu,2),size(obj.mu,3),class(obj.mu));
                    fprintf('Covariance \t(.Sigma) = \t[%i x %i x %i %s]\n',size(obj.Sigma,1),size(obj.Sigma,2),size(obj.Sigma,3),class(obj.Sigma))   ;
                    fprintf('Step Indexes(.kIndex)= \t[%i %s]\n',size(obj.kIndex,1),class(obj.kIndex))   ;
                end
                fprintf('Nr variables(.p) \t = \t%i\n',obj.p);
                fprintf('Nr steps (.l) \t \t = \t%i\n',obj.l);
                fprintf('Index method (.kMethod) = \t%s\n',obj.kMethod);
                disp(' ');
            else
                disp('Empty Time Variant Gaussian Noise Model.')
            end    
            
        end
        
        function [covmat,isDiag] = cov(obj,k)
        %COV Returns the covariance matrix.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_ltv object     
        %  K: step number
        %
        %  - Output variable(s) -
        %  COVMAT: covariance matrix Sigma at step K. If ISDIAG is not 
        %  supplied as output argument, the matrix is returned. If ISDIAG 
        %  is supplied then the original Sigma (matrix or vector) is returned.
        %
        %  ISDIAG: indicates if the covariance matrix is diagonal. 
        %  * ISDIAG = 1 : diagonal matrix, isDiag = 2 : diagonal vector, 
        %  * ISDIAG = 0 : none of the above         
        %
        %  - Construction -        
        %  [COVMAT]=COV(OBJ,K) returns the covariance matrix at step K.
        %
        %  [COVMAT,ISDIAG]=COV(OBJ,K) returns the covariance matrix or  
        %  vector at step K and ISDIAG.
        %
        %  [COVMAT]=COV(OBJ) returns the complete SIGMA array.
        %
        %  [COVMAT,ISDIAG]=COV(OBJ) returns the complete SIGMA array and 
        %  ISDIAG=0.        
        
            ni=nargin;
            
            if ni==1
                covmat=obj.Sigma;
                isDiag=0;
            elseif ni==2
                if size(obj.Sigma,3)>1      %Sigma is 3D array: find Sigma
                    covmat = findArrayVal(obj.Sigma,obj.kIndex,obj.kMethod,k);
                else                        %Sigma is matrix
                    covmat = obj.Sigma;
                end                
                
                if nargout>1                    %isDiag is needed
                    isDiag= isdiag(covmat); 
                else                            %if isDiag not needed: return full matrix
                    if min(size(covmat))==1 
                        covmat =diag(covmat);
                    end
                end
                
            else
                error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:argMismatch','Incorrect number of input arguments.');
            end

        end
 
        function [cholmat,isDiag] = cholcov(obj,k)
        %CHOLCOV Returns the upper triangular square root matrix of a 
        %        Positive-Semidefinite Sigma.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_ltv object     
        %
        %  - Output variable(s) -
        %  CHOLCOVMAT: upper triangular square root of covariance matrix 
        %  Sigma. If Sigma is semi-definite, the matrix is not triangular.  
        %  If ISDIAG is not supplied as output argument, the full matrix is 
        %  returned. If ISDIAG is supplied then the square root of the 
        %  original Sigma (matrix or vector) is returned.
        %
        %  ISDIAG: indicates if the covariance matrix is diagonal. 
        %  * ISDIAG = 1 : diagonal matrix, isDiag = 2 : diagonal vector, 
        %  * ISDIAG = 0 : none of the above   
        %
        %  - Construction -        
        %  [CHOLMAT]=COV(OBJ) returns the upper triangular square root
        %  matrix of Sigma.
        %
        %  [CHOLMAT]=COV(OBJ,~) returns the upper triangular square root
        %  matrix of Sigma.        
        %
        %  [CHOLMAT,ISDIAG]=COV(OBJ) returns the upper triangular square 
        %  root matrix or vector of Sigma and ISDIAG.
        %
        %  [CHOLMAT,ISDIAG]=COV(OBJ,~) returns the upper triangular square 
        %  root matrix or vector of Sigma and ISDIAG.       

            ni=nargin;
            
            if ni==2
                [tmpSigma,isDiag]=cov(obj,k); 	%find Sigma         
                [SigmaChol,code] = cholSigma(tmpSigma,1,isDiag); %Cholesky decomposition: semi definite is allowed
                if code == 3    % Indefinite
                    error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:sigmaIndef','Sigma must be positive semi-definite')
                end                

                if isDiag==2        %Sigma is diag. vector
                    cholmat = diag(SigmaChol);
                else                %Sigma is matrix
                    cholmat = SigmaChol;
                end

                if nargout>1     	%isDiag is needed
                    cholmat=SigmaChol;
                end                        
                
            else
                error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:argMismatch','Incorrect number of input arguments.');
            end        
        
        end                
        
        function varvec = var(obj,k)
        %VAR Returns the variance vector.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_ltv object     
        %  K: step number
        %
        %  - Output variable(s) -
        %  VARVEC: variance vector of Sigma at step K.  
        %
        %  - Construction -        
        %  VARVEC=VAR(OBJ,K) returns the variance vector of Sigma at step K.
        
            ni=nargin;
            
            if ni==2
                if size(obj.Sigma,3)>1      %Sigma is 3D array: find Sigma
                    varvec = findArrayVal(obj.Sigma,obj.kIndex,obj.kMethod,k);
                else                        %Sigma is matrix
                    varvec = obj.Sigma;
                end
            else
                error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:argMismatch','Incorrect number of input arguments.');
            end

            if size(varvec,2)~=1             %if Sigma is matrix, return the diagonal
                varvec=diag(varvec);
            end

        end               
        
        function meanvec = mean(obj,k)
        %MEAN Returns the mean vector.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_ltv object     
        %  K: step number
        %
        %  - Output variable(s) -
        %  MEANVEC: mean vector mu at step K.  
        %
        %  - Construction -        
        %  MEANVEC=MEAN(OBJ,K) returns the mean vector mu at step K. 
        %
        %  MEANVEC=MEAN(OBJ) returns the complete mean vector mu.         
        
            ni=nargin;
            
            if ni==1                    % arg k is not supplied
            
                meanvec=obj.mu;
            
            elseif ni==2                % arg k is supplied
                
                if size(obj.mu,3)>1
                    meanvec = findArrayVal(obj.mu,obj.kIndex,obj.kMethod,k);
                else
                    meanvec = obj.mu;
                end
                
            else
                
                error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:argMismatch','Incorrect number of input arguments.');
                
            end
            
        end
        
        function samples = sample(obj,n,k)
        %SAMPLE Draw sample(s) from multivariate normal distribution.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_ltv object   
        %  N: amount of required samples
        %  K: step number
        %
        %  - Output variable(s) -
        %  SAMPLES: samples from distribution at step K, where each sample 
        %  is a column vector.
        %  
        %  - Construction -          
        %  SAMPLES=SAMPLE(OBJ) returns one sample of the multivariate normal 
        %  distribution with mean mu and covariance Sigma at step K=kIndex(1). 
        %  
        %  SAMPLES=SAMPLE(OBJ,N) returns N samples of the multivariate normal 
        %  distribution with mean mu and covariance Sigma at step K=kIndex(1). 
        %
        %  SAMPLES=SAMPLE(OBJ,N,K) returns N samples of the multivariate normal 
        %  distribution with mean mu and covariance Sigma at step K. 
        
            ni=nargin;
            
            if ni==1; n = 1; k = obj.kIndex(1);
            elseif ni==2; k = obj.kIndex(1);
            elseif ni==3;
            else error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:argMismatch','Incorrect number of input arguments.');
            end
        
            tmpMu=mean(obj,k);              %find mu             
            [tmpSigma,isDiag]=cov(obj,k); 	%find Sigma 
                
            [SigmaChol,code] = cholSigma(tmpSigma,1,isDiag); %Cholesky decomposition: semi definite is allowed
            if code == 3    % Indefinite
                error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:sigmaIndef','Sigma must be positive semi-definite')
            end        
            
            if isDiag==2        %if diagonal: element-by-element binary vector multiplication
                samples = bsxfun(@times,randn(obj.p,n),SigmaChol) + repmat(tmpMu,1,n);
            elseif isDiag==1    %if diagonal matrix: element-by-element binary vector multiplication
                samples = bsxfun(@times,randn(obj.p,n),diag(SigmaChol)) + repmat(tmpMu,1,n);                
            else                %if not diagonal: matrix multiplication
                samples =  SigmaChol' * randn(obj.p,n) + repmat(tmpMu,1,n);
            end
        end
        
        function densEst = pdf(obj, x,k)
        %PDF Returns density from distribution.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_ltv object   
        %  X: matrix where each column represents a density evaluation point.
        %  K: step number        
        %
        %  - Output variable(s) -
        %  DENSEST: density of the multivariate normal distribution at 
        %  step K, evaluated at each column of X.  
        %  
        %  - Construction -          
        %  PDF(OBJ,X) returns the density of the multivariate normal 
        %  distribution with mean mu and covariance Sigma, evaluated
        %  at each column of X at step K=kIndex(1). 
        %
        %  PDF(OBJ,X,K) returns the density of the multivariate normal 
        %  distribution with mean mu and covariance Sigma, evaluated
        %  at each column of X at step K.         
        
            ni=nargin;
            
            if ni==2; k = obj.kIndex(1);
            elseif ni==3 ;
            else error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:argMismatch','Incorrect number of input arguments.');
            end
            
            [p, n] = size(x);
            if(p ~= obj.p)
                error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:xSize','X has wrong variable(row) size for density calculation.')
            end
            
            tmpMu=mean(obj,k);              %find mu   
            [tmpSigma,isDiag]=cov(obj,k); 	%find Sigma   
            
            [SigmaChol,code] = cholSigma(tmpSigma,0,isDiag); %Cholesky decomposition: semi definite is NOT allowed
            if code ~= 0        % Indefinite     
                error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:sigmaIndef','Sigma must be positive definite for pdf calculation')
            end   
            
            tmpPdfMax = 1 / sqrt((2 * pi)^obj.p * (detSigma(SigmaChol))^2); 

            x = x - repmat(tmpMu, 1, n);

            if isDiag == 2                                              %diagonal vector
                densEst = tmpPdfMax * exp(-0.5 * (x.^2)' * (1./tmpSigma) );  
          	elseif isDiag == 1                                          %diagonal matrix: use diag vector
            	densEst = tmpPdfMax * exp(-0.5 * (x.^2)' * (1./diag(tmpSigma)) );             
            else                                                        %symmetric matrix
                %idea: calculate inverse using known cholesky factorisation
                %Ax=b, find x with A=U'U => U'Ux=b   => 1)solve U'y=b   2)solve Ux=y                    
                if n < (obj.p/100)          %small amount of pdfs: for loop is faster
                    densEst=zeros(n,1);
                 	for i=1:n
                        y=SigmaChol'\x(:,i);
                        x2=SigmaChol\y;
                        densEst(i) = tmpPdfMax * exp(-0.5 * x(:,i)' * x2);             
                  	end
                else                        %larger amount of pdfs: matrix calculation is faster
                	y=SigmaChol'\x;
                  	x2=SigmaChol\y;      
                 	densEst = diag(tmpPdfMax * exp(-0.5 * x' * x2)); 
                end
            end
        end
    
        function value=isempty(obj)
        %ISEMPTY Check if object is empty.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_ltv object   
        %
        %  - Output variable(s) -
        %  VALUE: indicates if object is empty (=1) or not empty (=0).
        %  
        %  - Construction -          
        %  VALUE = ISEMPTY(OBJ) checks if object is empty or not. The
        %  object is considered empty when mu, Sigma, kIndex or kMethod
        %  are empty.
               
            if (~isempty(obj.mu)&& ~isempty(obj.Sigma)&& ~isempty(obj.kIndex)&& ~isempty(obj.kMethod))
                value=0;
            else
                value=1;
            end

        end       
        
    end %methods
    
    methods (Access = private) 

        function obj=checkConsist(obj)
        %CHECKCONSIST Checks if object is consistent.
        
            if ~isempty(obj)

                if size(obj.mu,1)~=size(obj.Sigma,1)
                	error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:sigmaMuDimMismatch','The dimensions of Sigma and mu must agree')
                end

                if (length(size(obj.mu))==3) && (length(size(obj.Sigma))==3) && (size(obj.mu,3)~=size(obj.Sigma,3))
                    error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:sigmaMuDimMismatch','The dimensions of mu and Sigma must agree')
                end
                
                if (length(size(obj.mu))==3) && (size(obj.mu,3)~=length(obj.kIndex))
                    error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:tIndexMuDimMismatch','The dimensions of mu and kIndex must agree')
                end
                
                if (length(size(obj.Sigma))==3) && (size(obj.Sigma,3)~=length(obj.kIndex))
                    error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:tIndexSigmaDimMismatch','The dimensions of Sigma and kIndex must agree')
                end
                
                if (length(size(obj.mu))~=3) && (length(size(obj.Sigma))~=3)
                    error('DA:NoiseModels:nm_gauss_ltv:nm_gauss_ltv:sigmaMuDimMismatch','Either mu or Sigma must be 3D array')
                end                
                
            end
        
        end
        
    end %methods
end %classdef