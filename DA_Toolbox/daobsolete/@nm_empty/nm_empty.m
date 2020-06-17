classdef nm_empty < nm_

    properties (Access = public)
        p;          %Amount of variables  
    end
    
    properties (Dependent = true, SetAccess = private)
        mu;
        Sigma;
    end
       
    methods 
    
        function obj = nm_empty(size)
            
            narginchk(1,1);

            obj.p = size;
            
        end        
        
        function obj = set.p(obj,value)
            
            value=checkArgs(value,'p');          %check p                    
            obj.p=value;          
            obj=checkConsist(obj);
            
        end                      

        function value = get.mu(obj)
            
            value=zeros(obj.p,1);
            
        end

        function value = get.Sigma(obj)
            
            value=zeros(obj.p);
            
        end
          
        function display(obj)

        	disp(' ');
            disp([inputname(1),' = '])
            disp(' ');            
            if ~isempty(obj)
                disp('Empty Noise Model.')
                disp('------------------')
                disp('Properties:');
                if obj.p < 10
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
                disp('Empty Empty-Noise Model.')
            end
            
        end
        
        function value = cov(obj,~)
        
            value=obj.Sigma;
            
        end
 
        function value = mean(obj,~)
        
            value=obj.mu;
            
        end
        
        function samples = sample(obj,n,~)
           
            if nargin==1; n = 1;
            elseif nargin==2;
            elseif nargin==3;    
            else error('DA:NoiseModels:nm_empty:nm_empty:argMismatch','Incorrect number of input arguments.');
            end
            
            samples = repmat(obj.mu,1, n);
        end
        
        function densEst = pdf(obj, x,~)

            if nargin<2;
            	error('DA:NoiseModels:nm_empty:nm_empty:argMismatch','Incorrect number of input arguments.');
            end
            
            [p, n] = size(x);

            if(p ~= obj.p)
                error('DA:NoiseModels:nm_empty:nm_empty:xSize','X has wrong variable(row) size for density calculation.')
            end

            densEst = repmat(obj.mu,1, n);
            
        end
        
        function value=isempty(obj)

            if (~isempty(obj.mu)&& ~isempty(obj.Sigma))
                value=0;
            else
                value=1;
            end

        end               

    end %methods
    
    methods (Access = private) 
    
        function obj=checkConsist(obj)
    
            if ~isempty(obj)

                if length(obj.mu)~=length(obj.Sigma)
                	error('DA:NoiseModels:nm_empty:nm_empty:sigmaMuDimMismatch','The dimensions of Sigma and mu must agree')
                end

            end
        
        end

    end %methods
end %classdef