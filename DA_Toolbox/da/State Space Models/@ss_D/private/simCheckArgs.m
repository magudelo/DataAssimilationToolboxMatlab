function arg=simCheckArgs(arg,strCode)
%SIMCHECKARGS Check simulation methods arguments 
%
%  - Input variable(s) -
%  ARG: argument to be checked
%  STRCODE: string code that specifies the type of argument.
%
%  - Output variable(s) -
%  ARG: checked and possibly adjusted argument
%
%  - Construction -
%  ARG = SIMCHECKARGS(ARG,STRCODE) checks the argument provided in ARG based
%  on the string code specified by STRCODE. If ARG does not contain a valid
%  format or value an error is generated.

    switch strCode
            
        case {'u'}              
            
            if ~isa(arg,'double')
                error('DA:StateSpaceModels:ss_D:checkArgs:ClassMismatch','%s must be of class ''double''.',strCode)
            end        

            N=size(arg);

            if length(N)>2                          %check dimensions
                error('DA:StateSpaceModels:ss_D:checkArgs:DimMismatch','%s must be a two-dimensional matrix.',strCode)
            end   

        case 'x0'
            
            x0=arg;
            
            if ~isa(x0,'double')
                error('DA:StateSpaceModels:ss_D:checkArgs:x0ClassMismatch','x0 must be of class ''double''.')
            end

            N=size(x0);

            if length(N)>2||(N(1)>1&&N(2)>1)        %check dimensions
                error('DA:StateSpaceModels:ss_D:checkArgs:x0DimMismatch','x0 must be a one dimensional column vector.')
            end

            if N(2)>1                               %check if column vector
                warning('DA:StateSpaceModels:ss_D:checkArgs:x0ColVec','x0 is transposed to a column vector.')
                x0=x0';
            end          
            
            arg=x0;                
              
        case 'samples'
            
            samples=arg;
            
            if isempty(samples)
                error('DA:StateSpaceModels:ss_D:checkArgs:samplesEmpty','%s must be non-empty.',strCode)
            end                  
            if ~isa(samples,'double')
                error('DA:StateSpaceModels:ss_D:checkArgs:samplesClassMismatch','%s must be of class ''double''.',strCode)
            end
            n=size(samples);
            if length(n)>2||(n(1)>1||n(2)>1)        %check dimensions
                error('DA:StateSpaceModels:ss_D:checkArgs:samplesScal','%s must be a scalar.',strCode)
            end        
            if samples<=0
                error('DA:StateSpaceModels:ss_D:checkArgs:samplesNeg','%s must be strictly positive.',strCode)
            end                  
            arg=samples;                 
            
              
        case 'noise'
            
            noise=arg;
            
            if isempty(noise)
                error('DA:StateSpaceModels:ss_D:checkArgs:noiseEmpty','%s must be non-empty.',strCode)
            end           
            
            if ~isa(noise,'double')
                error('DA:StateSpaceModels:ss_D:checkArgs:noiseClassMismatch','%s must be of class ''double''.',strCode)
            end
            
            n=size(noise);
            if length(n)>2||(n(1)>1||n(2)>1)        %check dimensions
                error('DA:StateSpaceModels:ss_D:checkArgs:noiseScal','%s must be a scalar.',strCode)
            end        
            if noise~=0 && noise~=1 && noise~=2 && noise~=3
                warning('DA:StateSpaceModels:ss_D:checkArgs:noiseNeg','%s must have value 0, 1, 2 or 3: %s set to default 3',strCode,strCode)
                noise=3;
            end                  
            arg=noise;                 
                        
            
        case {'conf'}              
            
            conf=arg;
            
            if ~isa(conf,'char')
                error('DA:StateSpaceModels:ss_D:checkArgs:ClassMismatch','%s must be of class ''char''.',strCode)
            end        
      
            arg=conf;
            
    end       

end