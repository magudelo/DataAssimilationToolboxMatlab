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
        if ~isa(mu,'function_handle') 
        	error('DA:NoiseModels:nm_gauss_handle:checkArgs:muClassMismatch','Mu must be of class ''function handle''.')
        end

        try
            tmpMu = mu(1,1);
        catch err
        	error('DA:NoiseModels:nm_gauss_handle:checkArgs:muErr','Error in function mu: %s',err.message)
        end          

        N=size(tmpMu);
            
        if length(N)>2||(N(1)>1&&N(2)>1)        %check dimensions
        	error('DA:NoiseModels:nm_gauss_handle:checkArgs:muDimMismatch','Mu must be a one dimensional column vector.')
        end
            
        if N(2)>1                               %check if column vector
            error('DA:NoiseModels:nm_gauss_handle:checkArgs:muColVec','Mu should be a column vector.')
        end        
            
        arg=mu;
        
    elseif strcmp(strCode,'Sigma')
        Sigma=arg;
        if ~isa(Sigma,'function_handle') 
        	error('DA:NoiseModels:nm_gauss_handle:checkArgs:sigmaClassMismatch','Sigma must be of class ''function handle''.')
        end
        
        try
            tmpSigma = Sigma(1,1);
        catch err
        	error('DA:NoiseModels:nm_gauss_handle:checkArgs:SigmaErr','Error in function Sigma: %s',err.message)
        end          

        N=size(tmpSigma);

        if length(N)>2                                  %check dimensions
        	error('DA:NoiseModels:nm_gauss_handle:checkArgs:sigmaDimMismatch','Sigma must be a two-dimensional matrix.')
        end
        
        if min(N)>1 && ~issym(tmpSigma)                 %if matrix, must be symmetric
            error('DA:NoiseModels:nm_gauss_handle:checkArgs:sigmaSym','Sigma must be square and symmetric')
        end
        
        if min(N)==1 && N(2)>1                        	%check if column vector
            error('DA:NoiseModels:nm_gauss_handle:checkArgs:sigmaColVec','Sigma must be a column vector.')
        end     
        
        arg=Sigma;

    elseif strcmp(strCode,'Ts')
            
            Ts=arg;
            
            if isempty(Ts)
                error('DA:NoiseModels:nm_gauss_handle:checkArgs:TsEmpty','%s must be non-empty.',strCode)
            end                  
            if ~isa(Ts,'double')
                error('DA:NoiseModels:nm_gauss_handle:checkArgs:TsClassMismatch','%s must be of class ''double''.',strCode)
            end
            n=size(Ts);
            if length(n)>2||(n(1)>1||n(2)>1)        %check dimensions
                error('DA:NoiseModels:nm_gauss_handle:checkArgs:TsScal','%s must be a scalar.',strCode)
            end        
            if Ts<=0
                error('DA:NoiseModels:nm_gauss_handle:checkArgs:TsNeg','%s must be strictly positive.',strCode)
            end                  
            arg=Ts;       
    
	end
end