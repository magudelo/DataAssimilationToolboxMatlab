function obj = nmAutoSel(varargin)
%NMAUTOSEL Automatic selection of noise model
%
%  - Input variable(s) -
%  TYPE: a string that defines the type of distribution. Possible values
%  are 'gauss' and 'empty'.
%
%  SUBTYPE: a cell or double that defines the subtype of the distribution.
%  For type 'gauss':
%  * double: when a 3D array is provided a LTV model is created, when a
%  regular array is provided a LTI model is created.
%  * cell: contains the input arguments of the noise model to create. For
%  example: {@funmu,@funSigma} creates nm_gauss_handle(@funmu,@funSigma).
%
%  - Output variable(s) -
%  OBJ: noise model object.
%
%  - Construction -
%  OBJ=NMAUTOSEL(TYPE,SUBTYPE) returns a noise model object as specified in
%  TYPE and SUBTYPE.
%
%  OBJ=NMAUTOSEL(SUBTYPE) returns a gaussian noise model object with
%  subtype SUBTYPE.
%
%  OBJ=NMAUTOSEL() returns the univariate standard normal distribution.
%  ( nm_gauss_lti( ) )

    subType=0;
	narginchk(0,2);   
	ni=nargin;
            
	switch ni
        case 0      %no arguments: LTI univariate standard normal distribution
            args={};            
            type='gauss';
        case 1      %one argument: subtype defined by args, type: normal distribution
            args=varargin{1};
            type='gauss';
                    
        case 2      %two arguments: args + type from arguments                    
            args=varargin{1};
            type=varargin{2};
	end
                        
	if isa(args,'cell')     %if cell: check content to define subtype
        
        for i=1:length(args)
                       
            if isa(args{i},'function_handle')   %nm_xxx_handle
                subType=3;
                break;
            end
            if isa(args{i},'double')            
                if length(size(args{i}))>2   	%nm_xxx_LTV
                	subType=2;
                    break;
                else                            %nm_xxx_LTI
                	subType=1;
                end
            end              
            
        end %for
        
        if isempty(args)      %empty cell: LTI univariate standard normal distribution
            subType=4;
        end
        
    %note: function handle is only possible in cell because two handles are
    %required.
        
    elseif isa(args,'double')       
        if length(size(args))>2             %nm_xxx_LTV
        	subType=2;
        else                                %nm_xxx_LTI
        	subType=1;
        end

	end %if
            
	if strcmp(type,'gauss')             %nm_gauss_xxx
        
        switch subType
            case 1
               	obj=nm_gauss_lti(args);
            case 2
               	obj=nm_gauss_ltv(args);
            case 3
                obj=nm_gauss_handle(args);
            case 4
                obj=nm_gauss_lti();                
            otherwise
                error('DA:NoiseModels:nm:argsErr','Invalid nm arguments.')
        end
        
    else
        error('DA:NoiseModels:nm:typeErr','Distribution type unknown.')
	end
            
end %nmAutoSel    