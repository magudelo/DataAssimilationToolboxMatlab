classdef nm_gauss_lti < nm_gauss
    %NM_GAUSS_LTI Constructs a multivariate Linear Time Invariant Gaussian 
    %Noise Model object
    %
    %  NM_GAUSS_LTI Properties:
    %     mu - Mean vector (Column vector)
    %     Sigma - Covariance matrix (Symmetric matrix or column vector that 
    %     represents the diagonal elements of a diagonal matrix.)
    %     p - amount of variables (read only)    
    %
    %  NM_GAUSS_LTI Construction:
    %     OBJ = NM_GAUSS_LTI(MU,SIGMA) creates an object OBJ representing
    %     a multivariate linear time invariant gaussian noise model with  
    %     mean MU and covariance matrix SIGMA.
    %
    %     OBJ = NM_GAUSS_LTI(SIGMA) creates a multivariate LTI gaussian 
    %     noise model with zero mean and covariance matrix SIGMA. 
    %     SIGMA can not be a diagonal vector in this case.
    %
    %     OBJ = NM_GAUSS_LTI(MU) creates a multivariate LTI gaussian noise
    %     model with mean MU and the identity matrix as covariance matrix.
    %     
    %     OBJ = NM_GAUSS_LTI() creates a univariate LTI gaussian noise model
    %     with zero mean and variance one (the standard normal distribution).
    %     
    %     OBJ = NM_GAUSS_LTI(CELL), with CELL a cell array, creates a 
    %     LTI gaussian noise model based on the information in CELL.
    %     
    %  NM_GAUSS_LTI Methods:
    %     cov - Returns the covariance matrix.
    %     var - Returns the variance matrix.
    %     mean - Returns the mean vector.
    %     sample - Draw sample(s) from distribution.
    %     pdf - Returns density from distribution.
    %     cholcov - Returns the upper triangular square root matrix of a 
    %               Positive-Semidefinite Sigma.    
    %
    %  See also NM_GAUSS_LTV, NM_GAUSS_HANDLE
    
    properties (Access = public)
        mu;         % Mean vector
        Sigma;      % Covariance matrix
    end
    
    properties (Dependent = true, SetAccess = private)
        p;          % Amount of variables    
    end
    
    properties (Access = private)
        SigmaChol;          %Cholesky factor of Sigma
        SigmaSemi;          %1= Is semi-positive definite
        pdfMax;             %precalculated constant
        isDiag;             %1= diagonal matrix, 2= diagonal vector, 0= neither
    end
    
    methods 
    
        function obj = nm_gauss_lti(varargin)
        % NM_GAUSS_LTI Constructor   
            narginchk(0, 2);
            ni = nargin;

            if ni==1 && isa(varargin{1},'cell')
                varargin=varargin{1};
                if length(size(varargin))> 2 || (size(varargin,1)>1 && size(varargin,2)>1) 
                    error('DA:NoiseModels:nm_gauss_lti:cellDim','Cell must be vector sized.')
                else
                    ni=max(size(varargin,1),size(varargin,2));
                end
            end
            
            switch ni
                case 0          
                    obj.mu=0;
                    obj.Sigma=1;
                case 1                   
                    if size(varargin{1},1) == size(varargin{1},2)   % if square matrix: Sigma
                        obj.Sigma=varargin{1};
                        obj.mu=zeros(obj.p,1);
                    else                                            % else: mu
                        obj.mu=varargin{1};
                        obj.Sigma=ones(obj.p,1);
                    end
                case 2                                              
                    obj.mu=varargin{1};                     
                    obj.Sigma=varargin{2};
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
            obj=calcConsts(obj);

        end                   
        
        function value = get.p(obj)                 %amount of variables 
            
            if isempty(obj.Sigma)
                value = length(obj.mu);
            else
                value = length(obj.Sigma);
            end
            
        end
          
        function display(obj)
        %DISPLAY Displays the object summary.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_lti object     
        %
        %  - Construction -
        %  DISPLAY(OBJ) displays the nm_gauss_lti object summary.
        
            if ~isempty(obj)              
                disp('Gaussian Noise Model - Time Invariant.')
                disp('--------------------------------------')
                disp('Properties:');
                if obj.p < 5
                    fprintf('Mean \t\t(.mu) \t = \n');
                    disp(obj.mu);
                    fprintf('Covariance \t(.Sigma) = \n');
                    disp(obj.Sigma);
                else
                    fprintf('Mean \t\t(.mu) \t = \t[%i x %i %s] \n',size(obj.mu,1),size(obj.mu,2),class(obj.mu))
                    fprintf('Covariance \t(.Sigma) = \t[%i x %i %s]\n',size(obj.Sigma,1),size(obj.Sigma,2),class(obj.Sigma))     
                end
                fprintf('Nr variables(.p) \t = \t%i\n',obj.p);
                disp(' ');
            else
                disp('Empty Time Invariant Gaussian Noise Model.')
            end                
        end
        
        function [covmat,isDiag] = cov(obj,~)
        %COV Returns the covariance matrix.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_lti object     
        %
        %  - Output variable(s) -
        %  COVMAT: covariance matrix Sigma. If ISDIAG is not supplied as 
        %  output argument, the matrix is returned. If ISDIAG is supplied
        %  then the original Sigma (matrix or vector) is returned.
        %
        %  ISDIAG: indicates if the covariance matrix is diagonal. 
        %  * ISDIAG = 1 : diagonal matrix, isDiag = 2 : diagonal vector, 
        %  * ISDIAG = 0 : none of the above         
        %
        %  - Construction -        
        %  [COVMAT]=COV(OBJ) returns the covariance matrix Sigma.
        %
        %  [COVMAT,ISDIAG]=COV(OBJ) returns the covariance matrix or vector
        %  and ISDIAG.
            
            isDiag=obj.isDiag;  %pre-calculated so no effort
            
            if isDiag==2        %Sigma is diag. vector
                covmat=diag(obj.Sigma);
            else                %Sigma is matrix
                covmat = obj.Sigma;
            end
                
            if nargout>1
                covmat=obj.Sigma;
            end
            
        end
 
        function [cholmat,isDiag] = cholcov(obj,~)
        %CHOLCOV Returns the upper triangular square root matrix of a 
        %        Positive-Semidefinite Sigma.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_lti object     
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

            isDiag=obj.isDiag;  %pre-calculated so no effort
            
            if isDiag==2        %Sigma is diag. vector
                cholmat = diag(obj.SigmaChol);
            else                %Sigma is matrix
                cholmat = obj.SigmaChol;
            end
                
            if nargout>1     	%isDiag is needed
                cholmat=obj.Sigma;
            end        
            
        end        
        
        function varvec = var(obj,~)
        %VAR Returns the variance vector.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_lti object     
        %
        %  - Output variable(s) -
        %  VARVEC: variance vector of Sigma.  
        %
        %  - Construction -        
        %  VARVEC=VAR(OBJ) returns the variance vector of Sigma.

            if obj.isDiag==2                %Sigma is diag. vector: var=vector
                varvec=obj.Sigma;
            else                            %Sigma is matrix: var= diagonal vector
                varvec=diag(obj.Sigma);
            end
            
        end        
        
        function meanvec = mean(obj,~)
        %MEAN Returns the mean vector.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_lti object     
        %
        %  - Output variable(s) -
        %  MEANVEC: mean vector mu.  
        %
        %  - Construction -        
        %  MEANVEC=MEAN(OBJ) returns the mean vector mu.             
        
            meanvec=obj.mu;
            
        end
        
        function samples = sample(obj,n,~)
        %SAMPLE Draw sample(s) from multivariate normal distribution.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_lti object   
        %  N: amount of required samples
        %
        %  - Output variable(s) -
        %  SAMPLES: samples from distribution, where each sample is a 
        %  column vector.
        %  
        %  - Construction -          
        %  SAMPLES=SAMPLE(OBJ) returns one sample of the multivariate normal 
        %  distribution with mean mu and covariance Sigma. 
        %  
        %  SAMPLES=SAMPLE(OBJ,N) returns N samples of the multivariate normal 
        %  distribution with mean mu and covariance Sigma. 
        
            if nargin==1; n = 1;
            elseif nargin==2;
            elseif nargin==3;    
            else error('DA:NoiseModels:nm_gauss_lti:nm_gauss_lti:argMismatch','Incorrect number of input arguments.');
            end

            if obj.isDiag==2        %if diagonal: element-by-element binary vector multiplication
                samples = bsxfun(@times,randn(obj.p,n),obj.SigmaChol) + repmat(obj.mu,1,n);
            elseif obj.isDiag==1    %if diagonal matrix: element-by-element binary vector multiplication
                samples = bsxfun(@times,randn(obj.p,n),diag(obj.SigmaChol)) + repmat(obj.mu,1,n);                
            else                    %if not diagonal: matrix multiplication
                samples =  obj.SigmaChol' * randn(obj.p,n) + repmat(obj.mu,1,n);
            end
            
        end
        
        function densEst = pdf(obj,x,~)
        %PDF Returns density from distribution.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_lti object   
        %  X: matrix where each column represents a density evaluation point.
        %
        %  - Output variable(s) -
        %  DENSEST: density of the multivariate normal distribution, evaluated
        %  at each column of X.  
        %  
        %  - Construction -          
        %  PDF(OBJ,X) returns the density of the multivariate normal 
        %  distribution with mean mu and covariance Sigma, evaluated
        %  at each column of X. 
        
            if nargin<2;
            	error('DA:NoiseModels:nm_gauss_lti:nm_gauss_lti:argMismatch','Incorrect number of input arguments.');
            end
            
            [p, n] = size(x);

            if(p ~= obj.p)
                error('DA:NoiseModels:nm_gauss_lti:nm_gauss_lti:xSize','X has wrong variable(row) size for density calculation.')
            end     
            
            x = x - repmat(obj.mu, 1, n);

            if obj.SigmaSemi
                error('DA:NoiseModels:nm_gauss_lti:nm_gauss_lti:sigmaSemi','Sigma must be positive definite for pdf calculation.')
            else
                if obj.isDiag == 2                                          %diagonal vector
                    densEst = obj.pdfMax * exp(-0.5 * (x.^2)' * (1./obj.Sigma) );  
                elseif obj.isDiag == 1                                      %diagonal matrix: use diag vector
                    densEst = obj.pdfMax * exp(-0.5 * (x.^2)' * (1./diag(obj.Sigma )) );             
                else                                                        %symmetric matrix
                    %idea: calculate inverse using known cholesky factorisation
                    %Ax=b, find x with A=U'U => U'Ux=b   => 1)solve U'y=b   2)solve Ux=y
                    if n < (obj.p/100)          %small amount of pdfs: for loop is faster
                        densEst=zeros(n,1);
                        for i=1:n
                            y=obj.SigmaChol'\x(:,i);
                            x2=obj.SigmaChol\y;
                            densEst(i) = obj.pdfMax * exp(-0.5 * x(:,i)' * x2);             
                        end
                    else                        %larger amount of pdfs: matrix calculation is faster
                        y=obj.SigmaChol'\x;
                        x2=obj.SigmaChol\y;      
                        densEst = diag(obj.pdfMax * exp(-0.5 * x' * x2)); 
                    end
                end
            end
        end
        
        function value=isempty(obj)
        %ISEMPTY Check if object is empty.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_lti object   
        %
        %  - Output variable(s) -
        %  VALUE: indicates if object is empty (=1) or not empty (=0).
        %  
        %  - Construction -          
        %  VALUE = ISEMPTY(OBJ) checks if object is empty or not. The
        %  object is considered empty when mu or Sigma are empty.
        %   
        
            if (~isempty(obj.mu)&& ~isempty(obj.Sigma))
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

                if length(obj.mu)~=length(obj.Sigma)
                	error('DA:NoiseModels:nm_gauss_lti:nm_gauss_lti:sigmaMuDimMismatch','The dimensions of Sigma and mu must agree')
                end

            end
        
        end
            
        function obj=calcConsts(obj)
        %CALCCONSTS Precalculate some intensive constants.
        
            obj.isDiag= isdiag(obj.Sigma); 
                
            [obj.SigmaChol,code] = cholSigma(obj.Sigma,1,obj.isDiag);
            if code == 0        % Positive definite
                obj.SigmaSemi = 0; 
            elseif code == 1    % Positive semi-definite
                obj.SigmaSemi = 1;
            else                % Indefinite
                error('DA:NoiseModels:nm_gauss_lti:nm_gauss_lti:sigmaIndef','Sigma must be positive semi-definite or positive definite')
            end
            obj.pdfMax = 1 / sqrt((2 * pi)^obj.p * (detSigma(obj.SigmaChol))^2); 
            
        end
        
    end %methods
end %classdef