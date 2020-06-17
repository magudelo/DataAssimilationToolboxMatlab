classdef nm_gauss_handle < nm_gauss
    %NM_GAUSS_HANDLE Constructs a multivariate Gaussian Noise Model object
    %
    %  NM_GAUSS_HANDLE Properties:
    %     mu - function handle to function with inputs step nr k and sample 
    %     time Ts and with output the mean column vector. (e.g. @funmean)
    %     Sigma - function handle to function with inputs step nr k and sample
    %     time Ts and with output the covariance matrix (or column vector in 
    %     case of uncorrelated matrix). (e.g. @funcovar) 
    %     Ts - sample time (default is 1)
    %     p - amount of variables (read only)    
    %
    %  NM_GAUSS_HANDLE Construction:
    %     OBJ = NM_GAUSS_HANDLE(MU,SIGMA,Ts) creates an object OBJ representing
    %     a multivariate gaussian noise model with sample time Ts.
    %
    %     OBJ = NM_GAUSS_HANDLE(MU,SIGMA) creates an object OBJ representing
    %     a multivariate gaussian noise model with sample time Ts=1.     
    %     
    %     OBJ = NM_GAUSS_HANDLE(CELL), with CELL a cell array, creates a 
    %     HANDLE gaussian noise model based on the information in CELL.
    %     
    %  NM_GAUSS_HANDLE Methods:
    %     cov - Returns the covariance matrix.
    %     var - Returns the variance matrix.
    %     mean - Returns the mean vector.
    %     sample - Draw sample(s) from distribution.
    %     pdf - Returns density from distribution.
    %     cholcov - Returns the upper triangular square root matrix of a 
    %               Positive-Semidefinite Sigma.      
    %
    %  See also NM_GAUSS_LTI, NM_GAUSS_LTV
    
    properties (Access = public)
        mu;             % function handle to obtain mean vector
        Sigma;          % function handle to obtain covariance matrix
        Ts;             % sample time (default=1)
    end
 
    properties (Dependent = true, SetAccess = private)
        p;              %Amount of variables    
    end    
    
    methods 
    
        function obj = nm_gauss_handle(varargin)
        % NM_GAUSS_HANDLE Constructor      
            ni=nargin;
            
            if ni==1 && isa(varargin{1},'cell')
                varargin=varargin{1};
                if length(size(varargin))> 2 || (size(varargin,1)>1 && size(varargin,2)>1) 
                    error('DA:NoiseModels:nm_gauss_handle:cellDim','Cell must be vector sized.')
                end
            else                
                narginchk(2, 3);    
            end
            
            obj.mu=varargin{1};
            obj.Sigma=varargin{2};               
            if length(varargin)<3
                obj.Ts=1;
            else
                obj.Ts=varargin{3};
            end
            if isempty(obj.Ts);obj.Ts=1;end;
            
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

        function obj = set.Ts(obj,value)
                                    
            value=checkArgs(value,'Ts');            %check Ts                    
            obj.Ts=value;
            obj=checkConsist(obj);
            
        end                 
        
        function value = get.p(obj)                 %amount of variables 
            tmpSigma = obj.Sigma(1,1);
            value = length(tmpSigma);
        end
        
        function display(obj)
        %DISPLAY Displays the object summary.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_handle object     
        %
        %  - Construction -
        %  DISPLAY(OBJ) displays the nm_gauss_handle object summary.
        
        	disp(' ');
            disp([inputname(1),' = '])
            disp(' ');            
            if ~isempty(obj)            
                disp('Gaussian Noise Model - Function Handle.')
                disp('---------------------------------------')
                disp('Properties:');
                fprintf('Mean \t\t(.mu) \t = %s\n',f2str(obj.mu));
                fprintf('Covariance \t(.Sigma) = %s\n',f2str(obj.Sigma));
                fprintf('Sample Time (.Ts) \t = %i \n',obj.Ts);                
                disp(' ');
            else
                disp('Empty Gaussian Noise - Function Handle Model.')
            end                
            
        end
        
        function [covmat,isDiag] = cov(obj,k)
        %COV Returns the covariance matrix.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_handle object     
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
        %  [COVMAT]=COV(OBJ,K) returns the covariance matrix at step number K.
        %
        %  [COVMAT,ISDIAG]=COV(OBJ,K) returns the covariance matrix or 
        %  vector at step number K and ISDIAG.
        %
        %  [COVMAT]=COV(OBJ) returns the function handle Sigma.
        %
        %  [COVMAT,ISDIAG]=COV(OBJ) returns the function handle Sigma and 
        %  ISDIAG=0. 
        
            ni=nargin;
            
            if ni==1
                covmat=obj.Sigma;
                isDiag = 0;
            elseif ni==2
                covmat=obj.Sigma(k,obj.Ts);
                
                if nargout>1                    %isDiag is needed
                    isDiag= isdiag(covmat); 
                else                            %if isDiag not needed: return full matrix
                    if min(size(covmat))==1 
                        covmat =diag(covmat);
                    end
                end                
                
            else
                error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:argMismatch','Incorrect number of input arguments.');
            end
            
        end
        
        function [cholmat,isDiag] = cholcov(obj,k)
        %CHOLCOV Returns the upper triangular square root matrix of a 
        %        Positive-Semidefinite Sigma.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_handle object     
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
                    error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:sigmaIndef','Sigma must be positive semi-definite')
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
                error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:argMismatch','Incorrect number of input arguments.');
            end        
        
        end                               
 
        function varvec = var(obj,k)
        %VAR Returns the variance vector.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_handle object     
        %  K: step number
        %
        %  - Output variable(s) -
        %  VARVEC: variance vector of Sigma step number K.  
        %
        %  - Construction -        
        %  VARVEC=VAR(OBJ,K) returns the variance vector of Sigma at step number K.
        
            ni=nargin;
            
            if ni==2
                varvec=obj.Sigma(k,obj.Ts);
            else
                error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:argMismatch','Incorrect number of input arguments.');
            end

            if size(varvec,2)~=1                %Sigma is matrix: var= diagonal vector
                varvec=diag(varvec);
            end

        end                
        
        function value = mean(obj,k)
        %MEAN Returns the mean vector.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_handle object     
        %  K: step number
        %
        %  - Output variable(s) -
        %  MEANVEC: mean vector mu at step number K.  
        %
        %  - Construction -        
        %  MEANVEC=MEAN(OBJ,K) returns the mean vector mu at step number K. 
        %
        %  MEANVEC=MEAN(OBJ) returns the function handle mu.     
        
            ni=nargin;
            
            if ni==1
            
                value=obj.mu;
            
            elseif ni==2
                
                value=obj.mu(k,obj.Ts);
                
            else
                
                error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:argMismatch','Incorrect number of input arguments.');
                
            end
            
        end
        
        function samples = sample(obj,n,k)
        %SAMPLE Draw sample(s) from multivariate normal distribution.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_handle object   
        %  N: amount of required samples
        %  K: step number
        %
        %  - Output variable(s) -
        %  SAMPLES: samples from distribution step number K, where each sample 
        %  is a column vector.
        %  
        %  - Construction -          
        %  SAMPLES=SAMPLE(OBJ,K) returns 1 sample of the multivariate normal 
        %  distribution with mean vector acquired from function handle mu 
        %  and covariance matrix acquired from function handle Sigma at  
        %  step number K.
        %
        %  SAMPLES=SAMPLE(OBJ,N,K) returns N samples of the multivariate  
        %  normal distribution with mean vector acquired from function handle  
        %  mu and covariance matrix acquired from function handle Sigma at  
        %  step number K=kIndex(1).
        
            if nargin==2; n = 1;
            elseif nargin==3;
            else error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:argMismatch','Incorrect number of input arguments.');
            end
            
            tmpMu=mean(obj,k);       %find mu     
            [tmpSigma,isDiag]=cov(obj,k); 	%find Sigma 
            
            [SigmaChol,code] = cholSigma(tmpSigma,1,isDiag); %Cholesky decomposition: semi definite is allowed
            if code == 3    % Indefinite
                error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:sigmaIndef','Sigma must be positive semi-definite')
            end        
            
           if isDiag==2        %if diagonal: element-by-element binary vector multiplication
                samples = bsxfun(@times,randn(obj.p,n),SigmaChol) + repmat(tmpMu,1,n);
            elseif isDiag==1    %if diagonal matrix: element-by-element binary vector multiplication
                samples = bsxfun(@times,randn(obj.p,n),diag(SigmaChol)) + repmat(tmpMu,1,n);                     
            else             	%if not diagonal: matrix multiplication
                samples =  SigmaChol' * randn(obj.p,n) + repmat(tmpMu,1,n);
            end
            
        end
        
        function densEst = pdf(obj,x,k)
        %PDF Returns density from distribution.
        %
        %  - Input variable(s) -
        %  OBJ: a nm_gauss_handle object   
        %  X: matrix where each column represents a density evaluation point.
        %  K: step number        
        %
        %  - Output variable(s) -
        %  DENSEST: density of the multivariate normal distribution at 
        %  step number K, evaluated at each column of X.  
        %  
        %  - Construction -          
        %  PDF(OBJ,X,K) returns the density of the multivariate normal 
        %  distribution with mean vector acquired from function handle mu 
        %  and covariance matrix acquired from function handle Sigma at  
        %  step number K and evaluated at each column of X. 
        
            if nargin~=3;
                error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:argMismatch','Incorrect number of input arguments.');
            end
                        
            [p, n] = size(x);
            if(p ~= obj.p)
                error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:xSize','X has wrong variable(row) size for density calculation.')
            end

            tmpMu=mean(obj,k);          %find mu 
            [tmpSigma,isDiag]=cov(obj,k); 	%find Sigma   
            
            [SigmaChol,code] = cholSigma(tmpSigma,0,isDiag); %Cholesky decomposition: semi definite is NOT allowed
            if code ~= 0        % Indefinite   
                error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:sigmaIndef','Sigma must be positive semi-definite')
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
        %  OBJ: a nm_gauss_handle object   
        %
        %  - Output variable(s) -
        %  VALUE: indicates if object is empty (=1) or not empty (=0).
        %  
        %  - Construction -          
        %  VALUE = ISEMPTY(OBJ) checks if object is empty or not. The
        %  object is considered empty when mu, Sigma or Ts are empty.
        %  
        
            if (~isempty(obj.mu)&& ~isempty(obj.Sigma)&& ~isempty(obj.Ts))
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

                tmpMu = obj.mu(1,obj.Ts);
                tmpSigma = obj.Sigma(1,obj.Ts);

                if length(tmpMu)~=length(tmpSigma)
                	error('DA:NoiseModels:nm_gauss_handle:nm_gauss_handle:sigmaMuDimMismatch','The dimensions of Sigma and mu must agree')
                end     
                
            end
        
        end
     
    end %methods
end %classdef