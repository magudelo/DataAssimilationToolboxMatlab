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

    switch strCode
    
        case {'f','h'}
           
            if isempty(arg)
                error('DA:SystemModels:ss_DNL:checkArgs:argEmpty','%s must be non-empty.',strCode)
            end
            
            if ~isa(arg,'function_handle')
                error('DA:SystemModels:ss_DNL:checkArgs:ClassMismatch','%s must be of class ''function_handle''.',strCode)
            end

        case {'fJacX','fJacW','hJacX','hJacV'}
               
            if ~isa(arg,'function_handle') && ~isempty(arg)
                error('DA:SystemModels:ss_DNL:checkArgs:ClassMismatch','%s must be of class ''function_handle''.',strCode)
            end            
            
        case 'x0'
            
            x0=arg;
            
            if isempty(x0)
                error('DA:SystemModels:ss_DNL:checkArgs:x0Empty','%s must be non-empty.',strCode)
            end         
            if ~isa(x0,'double')&&~isa(x0,'nm_')&&~isa(x0,'cell')
                error('DA:SystemModels:ss_DNL:checkArgs:x0ClassMismatch','Class mismatch %s.',strCode)
            end            
            arg=x0;                    

        case 'w'
            
            w=arg;

            if isempty(w)
                error('DA:SystemModels:ss_DNL:checkArgs:wEmpty','%s must be non-empty.',strCode)
            end             
            if ~isa(w,'double')&&~isa(w,'nm_')&&~isa(w,'cell')
                error('DA:SystemModels:ss_DNL:checkArgs:wClassMismatch','Class mismatch %s.',strCode)
            end            
            arg=w;                  
            
        case 'v'
            
            v=arg;
              
            if isempty(v)
                error('DA:SystemModels:ss_DNL:checkArgs:vEmpty','%s must be non-empty.',strCode)
            end                  
            if ~isa(v,'double')&&~isa(v,'nm_')&&~isa(v,'cell')
                error('DA:SystemModels:ss_DNL:checkArgs:vClassMismatch','Class mismatch %s.',strCode)
            end            
            arg=v;                              
            
        case 'Ts'
            
            Ts=arg;
            
            if isempty(Ts)
                error('DA:SystemModels:ss_DNL:checkArgs:TsEmpty','%s must be non-empty.',strCode)
            end                  
            if ~isa(Ts,'double')
                error('DA:SystemModels:ss_DNL:checkArgs:TsClassMismatch','%s must be of class ''double''.',strCode)
            end
            n=size(Ts);
            if length(n)>2||(n(1)>1||n(2)>1)        %check dimensions
                error('DA:SystemModels:ss_DNL:checkArgs:TsScal','%s must be a scalar.',strCode)
            end        
            if Ts<=0
                error('DA:SystemModels:ss_DNL:checkArgs:TsNeg','%s must be strictly positive.',strCode)
            end                  
            arg=Ts;          
            
        case 'k0'
            
            k0=arg;
            
            if isempty(k0)
                error('DA:SystemModels:ss_DNL:checkArgs:k0Empty','%s must be non-empty.',strCode)
            end                  
            if ~isa(k0,'double')
                error('DA:SystemModels:ss_DNL:checkArgs:k0ClassMismatch','%s must be of class ''double''.',strCode)
            end
            
            if length(k0)>1         %check dimensions
                error('DA:SystemModels:ss_DNL:checkArgs:k0Scal','%s must be a scalar.',strCode)
            end
            
            if rem(k0,1) ~= 0                       %remainder after division
                error('DA:SystemModels:ss_DNL:checkArgs:k0Int','%s must be an integer.',strCode)
            end            
            
            if k0<0
                error('DA:SystemModels:ss_DNL:checkArgs:k0Neg','%s must be non-negative.',strCode)
            end                  
            arg=k0;            
            
        case 'TimeUnit'
            
            TimeUnit=arg;
            
            if isempty(TimeUnit)
                error('DA:SystemModels:ss_DNL:checkArgs:TimeUnitEmpty','%s must be must be non-empty.',strCode)
            end
            if ~isa(TimeUnit,'char')
                error('DA:SystemModels:ss_DNL:checkArgs:TsClassMismatch','%s must be a string value.',strCode)
            end            
            
            if  ~strcmp(TimeUnit,'nanoseconds')&&~strcmp(TimeUnit,'microseconds')&&...
                ~strcmp(TimeUnit,'milliseconds')&&~strcmp(TimeUnit,'seconds')&&...
                ~strcmp(TimeUnit,'minutes')&&~strcmp(TimeUnit,'hours')&&...
                ~strcmp(TimeUnit,'days')&&~strcmp(TimeUnit,'weeks')&&...
                ~strcmp(TimeUnit,'months')&&~strcmp(TimeUnit,'years')
                
            warning('DA:SystemModels:ss_DNL:checkArgs:TimeUnitStringErr','%s contained wrong time unit. Set to default ''seconds''.',strCode)
            TimeUnit='seconds';
            
            end
            arg=TimeUnit;   
          
    end

end