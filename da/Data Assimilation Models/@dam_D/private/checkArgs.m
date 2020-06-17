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
        
        case 'alg'
            
            alg=arg;
            
            if isempty(alg)
                error('DA:DataAssimilationModels:dam_D:checkArgs:AlgEmpty','%s must be must be non-empty.',strCode)
            end
            if ~isa(alg,'char')
                error('DA:DataAssimilationModels:dam_D:checkArgs:AlgClassMismatch','%s must be a string value.',strCode)
            end            

            arg=alg;              
        
        case {'y','xf','u'}               
            
            if ~isa(arg,'double')
                error('DA:DataAssimilationModels:dam_D:checkArgs:ClassMismatch','%s must be of class ''double''.',strCode)
            end        

            N=size(arg);

            if length(N)>2                          %check dimensions
                error('DA:DataAssimilationModels:dam_D:checkArgs:DimMismatch','%s must be a two-dimensional matrix.',strCode)
            end                        

        case {'xa'}
    
            if  isempty(arg)
                error('DA:DataAssimilationModels:dam_D:checkArgs:xaEmpty','%s must be non-empty.',strCode)
            end                   
            
            if ~isa(arg,'double')
                error('DA:DataAssimilationModels:dam_D:checkArgs:xaClassMismatch','%s must be of class ''double''.',strCode)
            end        

            N=size(arg);

            if length(N)>2                          %check dimensions
                error('DA:DataAssimilationModels:dam_D:checkArgs:xaDimMismatch','%s must be a two-dimensional matrix.',strCode)
            end               

        case {'Pf','Pa','K'}
     
            if ~isa(arg,'double')
                error('DA:DataAssimilationModels:dam_D:checkArgs:ClassMismatch','%s must be of class ''double''.',strCode)
            end        

            N=size(arg);

            if length(N)>3                          %check dimensions
                error('DA:DataAssimilationModels:dam_D:checkArgs:DimMismatch','%s must be maximum a three-dimensional matrix.',strCode)
            end                
            
        case {'k'}
        
            k=arg;
            if  isempty(k)
                error('DA:SimulatorModels:dam_D:checkArgs:tEmpty','%s must be non-empty.',strCode)
            end                
            
            if ~isa(k,'double')
                error('DA:SimulatorModels:dam_D:checkArgs:tClassMismatch','%s must be of class ''double''.',strCode)
            end        

            N=size(k);

            if length(N)>2                          %check dimensions
                error('DA:SimulatorModels:dam_D:checkArgs:tDimMismatch','%s must be a vector.',strCode)
            end       
            
            if min(N(1),N(2))>1                     %check dimensions
                error('DA:SimulatorModels:dam_D:checkArgs:tDimMismatch','%s must be a vector.',strCode)
            end                
            
            if ~ismonotonic(k,1,'INCREASING')
                error('DA:SimulatorModels:dam_D:checkArgs:tIndexMono','%s must be strictly monotonically increasing.',strCode)
            end
            
            if ~isempty(find(k<0,1))
                error('DA:SimulatorModels:dam_D:checkArgs:tIndexNeg','%s can not contain negative time values.',strCode)
            end                  
            arg=k;     
            
        case 'Ts'
            
            Ts=arg;
            
            if isempty(Ts)
                error('DA:DataAssimilationModels:dam_D:checkArgs:TsEmpty','%s must be non-empty.',strCode)
            end                  
            if ~isa(Ts,'double')
                error('DA:DataAssimilationModels:dam_D:checkArgs:TsClassMismatch','%s must be of class ''double''.',strCode)
            end

            if length(Ts)>1
                error('DA:DataAssimilationModels:dam_D:checkArgs:TsScal','%s must be a scalar.',strCode)
            end        
            if Ts<=0
                error('DA:DataAssimilationModels:dam_D:checkArgs:TsNeg','%s must be strictly positive.',strCode)
            end                  
            arg=Ts;                
            
        case 'covFull'
            
            covFull=arg;
            
            if isempty(covFull)
                error('DA:DataAssimilationModels:dam_D:checkArgs:covFullEmpty','%s must be non-empty.',strCode)
            end          
            
            if ~isa(covFull,'double')
                error('DA:DataAssimilationModels:dam_D:checkArgs:covFullClassMismatch','%s must be of class ''double''.',strCode)
            end
            
            arg=covFull; 
            
        case 'TimeUnit'
            
            TimeUnit=arg;
            
            if isempty(TimeUnit)
                error('DA:DataAssimilationModels:dam_D:checkArgs:TimeUnitEmpty','%s must be must be non-empty.',strCode)
            end
            if ~isa(TimeUnit,'char')
                error('DA:DataAssimilationModels:dam_D:checkArgs:TsClassMismatch','%s must be a string value.',strCode)
            end            
            
            if  ~strcmp(TimeUnit,'nanoseconds')&&~strcmp(TimeUnit,'microseconds')&&...
                ~strcmp(TimeUnit,'milliseconds')&&~strcmp(TimeUnit,'seconds')&&...
                ~strcmp(TimeUnit,'minutes')&&~strcmp(TimeUnit,'hours')&&...
                ~strcmp(TimeUnit,'days')&&~strcmp(TimeUnit,'weeks')&&...
                ~strcmp(TimeUnit,'months')&&~strcmp(TimeUnit,'years')
                
                warning('DA:DataAssimilationModels:dam_D:checkArgs:TimeUnitStringErr','%s contained wrong time unit. Set to default ''seconds''.',strCode)
                TimeUnit='seconds';
            end
            arg=TimeUnit;               
            
    end       

end