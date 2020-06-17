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
        	error('DA:NoiseModels:nm_gauss_lti:checkArgs:muEmpty','Mu must be non-empty.')
        end
                
        if ~isa(mu,'double')
        	error('DA:NoiseModels:nm_gauss_lti:checkArgs:muClassMismatch','Mu must be of class ''double''.')
        end
        
        N=size(mu);
            
        if length(N)>2||(N(1)>1&&N(2)>1)        %check dimensions
        	error('DA:NoiseModels:nm_gauss_lti:checkArgs:muDimMismatch','Mu must be a one dimensional column vector.')
        end
            
        if N(2)>1                               %check if column vector
            warning('DA:NoiseModels:nm_gauss_lti:checkArgs:muColVec','Mu is transposed to a column vector.')
            mu=mu';
        end        
            
        arg=mu;
        
    elseif strcmp(strCode,'Sigma')
        Sigma=arg;
        if  isempty(Sigma)
        	error('DA:NoiseModels:nm_gauss_lti:checkArgs:sigmaEmpty','Sigma must be non-empty.')
        end
        
        if ~isa(Sigma,'double')
        	error('DA:NoiseModels:nm_gauss_lti:checkArgs:sigmaClassMismatch','Sigma must be of class ''double''.')
        end        
        
        N=size(Sigma);
        
        if length(N)>2                          %check dimensions
        	error('DA:NoiseModels:nm_gauss_lti:checkArgs:sigmaDimMismatch','Sigma must be a column vector or a two-dimensional matrix.')
        end
        
        if min(N)>1 && ~issym(Sigma)            %if matrix, must be symmetric
            error('DA:NoiseModels:nm_gauss_lti:checkArgs:sigmaSym','Sigma must be square and symmetric')
        end
        
        if min(N)==1 && N(2)>1                  %check if column vector
            warning('DA:NoiseModels:nm_gauss_lti:checkArgs:sigmaColVec','Sigma is transposed to a column vector.')
            Sigma=Sigma';
        end     
        
        arg=Sigma;        
    end

end