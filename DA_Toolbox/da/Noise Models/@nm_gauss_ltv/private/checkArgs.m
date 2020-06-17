function arg=checkArgs(arg,strCode)
%CHECKARGS Check arguments 
%
%  - Input variable(s) -
%  ARG: argument to be checked
%  STRCODE: string code that specifies the type of argument.
%
%  - Output variable(s) -
%  ARG: checked and possibly adjusted argument
%
%  - Construction -
%  ARG = CHECKARGS(ARG,STRCODE) checks the argument provided in ARG based
%  on the string code specified by STRCODE. If ARG does not contain a valid
%  format or value an error is generated.

    if strcmp(strCode,'mu')
        mu=arg;
        if  isempty(mu)
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:muEmpty','Mu must be non-empty.')
        end
                
        if ~isa(mu,'double')
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:muClassMismatch','Mu must be of class ''double''.')
        end
        
        N=size(mu);
            
        if length(N)>3||(N(1)>1&&N(2)>1)        %check dimensions
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:muDimMismatch','Mu must be a 3D column vector.')
        end
            
        if N(2)>1                               %check if column vector
            warning('DA:NoiseModels:nm_gauss_ltv:checkArgs:muColVec','Mu is transposed to a column vector.')
            mu=permute(mu,[2 1 3]);
        end        
            
        arg=mu;
        
    elseif strcmp(strCode,'Sigma')
        Sigma=arg;
        if  isempty(Sigma)
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:sigmaEmpty','Sigma must be non-empty.')
        end
        
        if ~isa(Sigma,'double')
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:sigmaClassMismatch','Sigma must be of class ''double''.')
        end        
        
        N=size(Sigma);
        
        if length(N)>3                          %check dimensions
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:sigmaDimMismatch','Sigma must be maximum a three-dimensional matrix.')
        end
        
        if min( [N(1) N(2)] ) > 1 && ~issym(Sigma)
            error('DA:NoiseModels:nm_gauss_ltv:checkArgs:sigmaSym','Sigma must be square and symmetric')
        end
        
        if min( [N(1) N(2)] ) == 1 && N(2)>1                	%check if column vector
            warning('DA:NoiseModels:nm_gauss_ltv:checkArgs:sigmaColVec','Sigma is transposed to a column vector.')
            Sigma=permute(Sigma,[2 1 3]);
        end     
        
        arg=Sigma;
        
    elseif strcmp(strCode,'kIndex')
        kIndex=arg;
        if  isempty(kIndex)
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:kIndexEmpty','kIndex must be non-empty.')
        end        
        if ~isa(kIndex,'double')
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:kIndexClassMismatch','kIndex must be of class ''double''.')
        end                
        
        N=size(kIndex);
        
        if length(N)>2||(N(1)>1&&N(2)>1)        %check dimensions
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:kIndexVec','kIndex must be a one dimensional vector.')
        end
        
        if ~ismonotonic(kIndex,1,'INCREASING')
            error('DA:NoiseModels:nm_gauss_ltv:checkArgs:kIndexMono','kIndex must be strictly monotonically increasing.')
        end
        
        if ~isempty(find(kIndex<0,1))
            error('DA:NoiseModels:nm_gauss_ltv:checkArgs:kIndexNeg','kIndex can not contain negative values.')
        end             
        
        arg=kIndex;
        
    elseif strcmp(strCode,'kMethod')
        kMethod=arg;
        if  isempty(kMethod)
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:kMethodEmpty','kMethod must be non-empty.')
        end            
        
        if ~isa(kMethod,'char')
        	error('DA:NoiseModels:nm_gauss_ltv:checkArgs:kMethodClassMismatch','kMethod must be of class ''char''.')
        end    
        
        if ~strcmp(kMethod,'low')&&~strcmp(kMethod,'high')&&~strcmp(kMethod,'near')
            warning('DA:NoiseModels:nm_gauss_ltv:checkArgs:kMethodStringMismatch','kMethod has wrong string value. Reversed to ''low''.')
            kMethod='low';
        end
        
        arg=kMethod;
    end
       

end