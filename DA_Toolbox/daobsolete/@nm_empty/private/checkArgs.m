function arg=checkArgs(arg,strCode)

    if strcmp(strCode,'p')
        p=arg;
        
        if  isempty(p)
        	error('DA:NoiseModels:nm_empty:checkArgs:muEmpty','p must be non-empty.')
        end
                
        if ~isa(p,'double')
        	error('DA:NoiseModels:nm_empty:checkArgs:muClassMismatch','p must be of class ''double''.')
        end
        
        N=size(p);
            
        if (N(1)>1||N(2)>1)        %check dimensions
        	error('DA:NoiseModels:nm_empty:checkArgs:muDimMismatch','p must be a scalar.')
        end
          
        arg=p;
   
    end

end