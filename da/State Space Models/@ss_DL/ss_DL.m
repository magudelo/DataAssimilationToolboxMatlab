classdef ss_DL < ss_D
    %SS_DL Constructs a Discrete Linear Time Variant State-Space object
    %
    %      x(k+1)  = A(k)x(k) + B(k)u(k) + w(k)
    %      y(k)    = C(k)x(k) + D(k)u(k) + v(k)    
    %
    %  SS_DL Properties:
    %     A - State matrix (xSize x xSize x l)
    %
    %     B - Input matrix (xSize x uSize x l)
    %    
    %     C - Output matrix (ySize x xSize x l)
    %
    %     D - Feedthrough matrix (ySize x uSize x l). (When not used can be 
    %         set equal to the scalar 0)
    %
    %     x0 - Initial state - noise model. Type 'help nm_' for all possible
    %     noise models. Other possibilities are to enter a cell array of 
    %     gaussian noise model input arguments or to enter an array which  
    %     automatically creates a corresponding gaussian noise model.
    %       Examples:
    %       - noise model: x0 = nm_gauss_handle(@fun_mu,@fun_sigma)
    %       - cell array:  x0 = {@fun_mu,@fun_sigma}
    %       - array: x0 = [0;0;0;0] => zero mean LTI with unit variance
    %                x0 = [2 0;0 2] => zero mean LTI with variance [2 0;0 2]
    %                x0 = 3D-array  => LTV model
    %
    %     w -  Process noise - noise model. See x0 for more info.
    %
    %     v -  Measurement noise - noise model. See x0 for more info.
    %
    %     k0 - initial step number. (default is 0)
    %
    %     Ts - sample time. (default is 1)
    %
    %     TimeUnit - String representing the unit of the time variable. 
    %       TimeUnit can take the following values: 'nanoseconds',
    %       'microseconds', 'milliseconds', 'seconds', 'minutes', 'hours',
    %       'days', 'weeks', 'months', 'years'.
    %       Default: 'seconds'
    %
    %     kIndex - step index vector used to indicate which step number
    %     corresponds to each vector or matrix in case of 3D arrays 
    %     (A,B,C and/or D). (Default = values '0' to 'l-1').
    %     For example:
    %     If A is the 2 x 2 x 3 array of [1 2;2 1], [2 3;3 2] and 
    %     [3 4;4 3] and kIndex=[2 3 5] then the matrix [2 3;3 2] is the
    %     the A matrix at k=3.
    %
    %     kMethod - string that indicates the rounding/interpolation method
    %     to retrieve the 2D array from the 3D array at a step k while taking
    %     into account kIndex. The following strings are possible (default='low'):
    %     1) 'low': find 2D array corresponding to the lower bound step nr in kIndex 
    %     2) 'high': find 2D array corresponding to the higher bound step nr in kIndex
    %     3) 'near': find 2D array corresponding to the nearest step nr in kIndex
    %
    %     uSize - amount of inputs (read only)  
    %
    %     xSize - amount of states (read only)  
    %
    %     ySize - amount of outputs (read only)  
    %
    %     l - length of 3D-array(s), the third dimension (read only)  
    %
    %  SS_DL Construction:
    %     SS_DLMODEL=SS_DL(A,B,C,D,X0,W,V,K0,TS,TIMEUNIT,KINDEX,KMETHOD)
    %     creates a Discrete Linear Time Variant State-Space object 
    %     SS_DLMODEL. When an empty object is inserted for K0, TS, TIMEUNIT 
    %     KINDEX or KMETHOD their default value is used.
    %
    %     SS_DLMODEL=SS_DL(A,B,C,D,X0,W,V,K0,TS,TIMEUNIT,KINDEX)
    %     uses KMETHOD='low'.
    %     
    %     SS_DLMODEL=SS_DL(A,B,C,D,X0,W,V,K0,TS,TIMEUNIT)
    %     uses values '0' to 'l-1' for KINDEX and KMETHOD='low'.
    %         
    %     SS_DLMODEL=SS_DL(A,B,C,D,X0,W,V,K0,TS) uses TIMEUNIT='seconds',
    %     values '0' to 'l-1' for KINDEX and KMETHOD='low'.
    %      
    %     SS_DLMODEL=SS_DL(A,B,C,D,X0,W,V,K0) uses TS=1, TIMEUNIT='seconds',
    %     values '0' to 'l-1' for KINDEX and KMETHOD='low'.
    %
    %     SS_DLMODEL=SS_DL(A,B,C,D,X0,W,V) uses K0=0, TIMEUNIT ='seconds',
    %     TS=1, values '0' to 'l-1' for KINDEX and KMETHOD='low'.
    %
    %  SS_DL Methods:
    %     eval_A - Retrieve value of matrix A.
    %     eval_B - Retrieve value of matrix B.
    %     eval_C - Retrieve value of matrix C.
    %     eval_D - Retrieve value of matrix D.
    %     eval_ftot - Calculates x(k+1) of state equation.
    %     eval_htot - Calculates y(k) of measurement equation.	
    %     eval_ftotJacX - Calculates Jacobian with respect to x of the 
    %                     state equation.
    %     eval_htotJacX - Calculates Jacobian with respect to x of the 
    %                     measurement equation.  
    %     eval_ftotJacW - Calculates Jacobian with respect to w of the 
    %                     state equation.
    %     eval_htotJacV - Calculates Jacobian with respect to v of the 
    %                     measurement equation.
    %
    %  See also SS_DNL, SS_DNL_AN
    
    properties (Access = public)

        A;          % State matrix (xSize x xSize x l)
        B;          % Input matrix (xSize x uSize x l)
        C;      	% Output matrix (ySize x xSize x l)
        D;       	% Feedthrough matrix (ySize x uSize x l)
        x0;         % Initial state - noise model.
        w;          % Process noise - noise model.
        v;          % Measurement noise - noise model.
        k0;         % initial step number. 
        Ts;         % sample time.
        TimeUnit;   % String representing the unit of the time variable. 
        kIndex;     % Step index array
        kMethod;    % String representing rounding/interpolation method
        
    end

    properties (Dependent = true, SetAccess = private)
        uSize;      % amount of inputs
        xSize;      % amount of states
        ySize;      % amount of outputs    
        l;          % length of 3D-array
    end        
    
    methods
        
        function obj = ss_DL(A,B,C,D,x0,w,v,k0,Ts,TimeUnit,kIndex,kMethod)
        % SS_DL Constructor    
        
            narginchk(7, 12);
            ni = nargin;

            if ni>=7
                    obj.A=A;
                    obj.B=B;
                    obj.C=C;                    
                    obj.D=D;
                    obj.x0=x0;
                    obj.w=w;
                    obj.v=v;
            end

            if ni<8;k0=0;end;            
            if ni<9;Ts=1;end;            
            if ni<10;TimeUnit='seconds';end;            
            if ni<11;kIndex=(0:obj.l-1);end;            
            if ni<12;kMethod='low';end;
            
            if isempty(k0);k0=0;end; 
            if isempty(Ts);Ts=1;end; 
            if isempty(TimeUnit);TimeUnit='seconds';end;            
            if isempty(kIndex);kIndex=(0:obj.l-1);end;            
            if isempty(kMethod);kMethod='low';end;

            obj.k0=k0;
            obj.Ts=Ts;
            obj.TimeUnit=TimeUnit;
            obj.kIndex=kIndex;
            obj.kMethod=kMethod;
            
        end
        
        function obj = set.A(obj,value)
            
            value=checkArgs(value,'A');             %check A   
            obj.A=value;
            checkConsist(obj);

        end         
        
        function obj = set.B(obj,value)
            
            value=checkArgs(value,'B');         	%check B             
            obj.B=value;
            checkConsist(obj);
            
        end         
        
        function obj = set.C(obj,value)
            
            value=checkArgs(value,'C');             %check C              
            obj.C=value;
            checkConsist(obj);
            
        end               
        
        function obj = set.D(obj,value)
            
            value=checkArgs(value,'D');             %check D   
            if isequal(value,0)
                value=zeros(obj.ySize,obj.uSize); %#ok<MCSUP>
            end
            obj.D=value;
            checkConsist(obj);
            
        end   
        
        function obj = set.x0(obj,value)
            
            value=checkArgs(value,'x0');            %check x0     
            if ~isa(value,'nm_')     
                obj.x0=nmAutoSel(value,'gauss');
            else
                obj.x0=value;
            end            
            checkConsist(obj);
            
        end          
        
        function obj = set.w(obj,value)
            
            value=checkArgs(value,'w');             %check w   
            if ~isa(value,'nm_')
                obj.w=nmAutoSel(value,'gauss');
            else
                obj.w=value;
            end
            checkConsist(obj);
            
        end   
        
        function obj = set.v(obj,value)
            
            value=checkArgs(value,'v');             %check v                    
            if ~isa(value,'nm_')
                obj.v=nmAutoSel(value,'gauss');
            else
                obj.v=value;
            end
            checkConsist(obj);
            
        end           
        
        function obj = set.k0(obj,value)
            
            value=checkArgs(value,'k0');            %check k0            
            obj.k0=value;
            checkConsist(obj);
            
        end            
        
        function obj = set.Ts(obj,value)
            
            value=checkArgs(value,'Ts');            %check Ts             
            obj.Ts=value;
            checkConsist(obj);
            
        end    
        
        function obj = set.TimeUnit(obj,value)
            
            value=checkArgs(value,'TimeUnit');      %check TimeUnit                   
            obj.TimeUnit=value;
            checkConsist(obj);
            
        end        
        
        function obj = set.kIndex(obj,value)
            
            value=checkArgs(value,'kIndex');        %check kIndex                   
            obj.kIndex=value;
            checkConsist(obj);
            
        end             
        
        function obj = set.kMethod(obj,value)
            
            value=checkArgs(value,'kMethod');       %check kMethod                   
            obj.kMethod=value;
            checkConsist(obj);
            
        end                          
        
        function value = get.xSize(obj)
            
            value = max([size(obj.A,1) size(obj.B,1) size(obj.C,2)]);
            
        end
        
        function value = get.uSize(obj)
            
            value = max([size(obj.B,2) size(obj.D,2)]);
            
        end        
        
        function value = get.ySize(obj)
            
            value = max([size(obj.C,1) size(obj.D,1)]);
            
        end             

        function value = get.l(obj)
            
            value = max([size(obj.A,3) size(obj.B,3) size(obj.C,3) size(obj.D,3)]);
            
        end        
        
        function display(obj)
        %DISPLAY Displays the object summary.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object     
        %
        %  - Construction -
        %  DISPLAY(OBJ) displays the ss_DL object summary.
        
        	disp(' ');
            disp([inputname(1),' = '])
            disp(' ');          
            if ~isempty(obj)              
                disp('Discrete Linear State-Space Model.')
                disp('----------------------------------')
                disp('Properties:');
                if obj.l==1
                    fprintf('A (.A) \t = \t[%i x %i %s] \n',size(obj.A,1),size(obj.A,2),class(obj.A));
                    fprintf('B (.B) \t = \t[%i x %i %s] \n',size(obj.B,1),size(obj.B,2),class(obj.B));
                    fprintf('C (.C) \t = \t[%i x %i %s] \n',size(obj.C,1),size(obj.C,2),class(obj.C));
                    fprintf('D (.D) \t = \t[%i x %i %s] \n',size(obj.D,1),size(obj.D,2),class(obj.D));
                else
                    fprintf('A (.A) \t = \t[%i x %i x %i %s] \n',size(obj.A,1),size(obj.A,2),size(obj.A,3),class(obj.A));
                    fprintf('B (.B) \t = \t[%i x %i x %i %s] \n',size(obj.B,1),size(obj.B,2),size(obj.B,3),class(obj.B));
                    fprintf('C (.C) \t = \t[%i x %i x %i %s] \n',size(obj.C,1),size(obj.C,2),size(obj.C,3),class(obj.C));
                    fprintf('D (.D) \t = \t[%i x %i x %i %s] \n',size(obj.D,1),size(obj.D,2),size(obj.D,3),class(obj.D));
                end
                fprintf('Initial State (.x0)\t = \t[%s] \n',class(obj.x0));
                fprintf('Process noise (.w) \t = \t[%s] \n',class(obj.w));
                fprintf('Measur. noise (.v) \t = \t[%s] \n',class(obj.v));
                fprintf('Nr. States (.xSize)  = \t%i \n',obj.xSize);
                fprintf('Nr. Outputs(.ySize)  = \t%i \n',obj.ySize);
                fprintf('Nr. Inputs (.uSize)  = \t%i \n',obj.uSize);
                fprintf('Initial step nr (.k0)= \t%i \n',obj.k0);                
                fprintf('Sample Time(.Ts) \t = \t%i \n',obj.Ts);
                fprintf('Time Unit (.TimeUnit)= \t%s \n',obj.TimeUnit);
                if obj.l ~= 1
                    fprintf('Step Indexes(.kIndex) = [%i x %i %s] \n',size(obj.kIndex,1),size(obj.kIndex,2),class(obj.kIndex));
                    fprintf('Index Method(.kMethod)= %s \n',obj.kMethod);
                end
                disp(' ');
            else
                disp('Empty Discrete Linear State-Space Model.')
            end                
        end        
        
        function value=isempty(obj)
        %ISEMPTY Check if object is empty.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object   
        %
        %  - Output variable(s) -
        %  VALUE: indicates if object is empty (=1) or not empty (=0).
        %  
        %  - Construction -          
        %  VALUE = ISEMPTY(OBJ) checks if object is empty or not. The
        %  object is considered empty when A, B, C, D, x0, w, v, Ts, 
        %  TimeUnit, kIndex or kMethod are empty.
        
            if (~isempty(obj.A)&& ~isempty(obj.B)&& ~isempty(obj.C)&& ~isempty(obj.D)&&...
                ~isempty(obj.x0)&& ~isempty(obj.w)&&~isempty(obj.v)&& ~isempty(obj.Ts)&& ...
                ~isempty(obj.TimeUnit)&& ~isempty(obj.kIndex)&& ~isempty(obj.kMethod))
                value=0;
            else
                value=1;
            end

        end                
        
        function value=eval_A(obj,k)
        %EVAL_A Retrieve value of matrix A.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object   
        %  K: step number
        %
        %  - Output variable(s) -
        %  VALUE: matrix A (at step K).
        %  
        %  - Construction -          
        %  VALUE = EVAL_A(OBJ,K) retrieves matrix A at step K.
        %
        %  VALUE = EVAL_A(OBJ) retrieves matrix A. Can not be used if 
        %  A is a 3D array.       
        
            if size(obj.A,3)>1
                value = findArrayVal(obj.A,obj.kIndex,obj.kMethod,k);
            else
                value = obj.A;  
            end
            
        end
    
        function value=eval_B(obj,k)
        %EVAL_B Retrieve value of matrix B.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object   
        %  K: step number
        %
        %  - Output variable(s) -
        %  VALUE: matrix B (at step K).
        %  
        %  - Construction -          
        %  VALUE = EVAL_B(OBJ,K) retrieves matrix B at step K.
        %
        %  VALUE = EVAL_B(OBJ) retrieves matrix B. Can not be used if 
        %  B is a 3D array.    
        
            if size(obj.B,3)>1
                value = findArrayVal(obj.B,obj.kIndex,obj.kMethod,k);
            else
                value = obj.B;  
            end
            
        end        
        
        function value=eval_C(obj,k)
        %EVAL_C Retrieve value of matrix C.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object   
        %  K: step number
        %
        %  - Output variable(s) -
        %  VALUE: matrix C (at step K).
        %  
        %  - Construction -          
        %  VALUE = EVAL_C(OBJ,K) retrieves matrix C at step K.
        %
        %  VALUE = EVAL_C(OBJ) retrieves matrix C. Can not be used if 
        %  C is a 3D array. 
        
            if size(obj.C,3)>1
                value = findArrayVal(obj.C,obj.kIndex,obj.kMethod,k);
            else
                value = obj.C;  
            end
            
        end             
     
        function value=eval_D(obj,k)
        %EVAL_D Retrieve value of matrix D.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object   
        %  K: step number
        %
        %  - Output variable(s) -
        %  VALUE: matrix D (at step K).
        %  
        %  - Construction -          
        %  VALUE = EVAL_D(OBJ,K) retrieves matrix D at step K.
        %
        %  VALUE = EVAL_D(OBJ) retrieves matrix D. Can not be used if 
        %  D is a 3D array.
        
            if size(obj.D,3)>1
                value = findArrayVal(obj.D,obj.kIndex,obj.kMethod,k);
            else
                value = obj.D;  
            end
            
        end           
        
        function value=eval_ftot(obj,x,k,u,w)
        %EVAL_FTOT Calculates x(k+1) of state equation.
        %   x(k+1) = A(k)x(k) + B(k)u(k) + w(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object 
        %  X: state at step K (X=[] or X=0 calculates B(k)u(k) + w(k))
        %  K: step number
        %  U: input at step K (U=[] or X=0 calculates A(k)x(k) + w(k))
        %  W: process noise at step K. (W=[] takes sample from distribution, 
        %     W=0 calculates A(k)x(k) + B(k)u(k) )
        %
        %  - Output variable(s) -
        %  VALUE: x(k+1) of state equation.
        %  
        %  - Construction -          
        %  VALUE = EVAL_FTOT(OBJ,X,K,U,W) calculates x(k+1).
        %
        %  VALUE = EVAL_FTOT(OBJ,X,K,U) calculates x(k+1) with process 
        %  noise W a sample from its distribution.
        
            narginchk(4, 5);
            ni = nargin;  
            
            if ni~=5;w=[];end;
            
            value=zeros(obj.xSize,1);
            
            if ~isequal(x,0) && ~isempty(x)
                value = value + eval_A(obj,k) * x;   
            end 
               
            if ~isequal(u,0) && ~isempty(u)
                value = value + eval_B(obj,k)*u;
            end
            
            if isempty(w)               %if w=[]: take sample
                w=sample(obj.w,1,k);
                value = value + w;
            elseif ~isequal(w,0)        %if w=0: add nothing
                value = value + w;                
            end
            	
        end          
        
        function value=eval_htot(obj,x,k,u,v)
        %EVAL_HTOT Calculates y(k) of measurement equation.
        %   y(k) = C(k)x(k) + D(k)u(k) + v(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object 
        %  X: state at step K (X=[] or X=0 calculates D(k)u(k) + v(k))
        %  K: step number
        %  U: input at step K (U=[] or X=0 calculates C(k)x(k) + v(k))
        %  V: measurement noise at step K. (V=[] takes sample from 
        %     distribution, V=0 calculates C(k)x(k) + D(k)u(k) )
        %
        %  - Output variable(s) -
        %  VALUE: y(k) of measurement equation.
        %  
        %  - Construction -          
        %  VALUE = EVAL_HTOT(OBJ,X,K,U,V) calculates y(k).
        %
        %  VALUE = EVAL_HTOT(OBJ,X,K,U) calculates y(k) with measurement 
        %  noise V a sample from its distribution.
        
            narginchk(4, 5);
            ni = nargin;  
            
            if ni~=5;v=[];end;
            
            value=zeros(obj.ySize,1);
            
            if ~isequal(x,0) && ~isempty(x) 
                value = value + eval_C(obj,k) * x;    
            end 
               
            if ~isequal(u,0) && ~isempty(u) 
                value = value + eval_D(obj,k)*u;    
            end
            
            if isempty(v)               %if v=[]: take sample
                v=sample(obj.v,1,k);
                value = value + v;
            elseif ~isequal(v,0)        %if v=0: add nothing
                value = value + v;                
            end
            	
        end           

        function value=eval_ftotJacX(obj,~,k,~,~)
        %EVAL_FTOTJACX Calculates Jacobian with respect to x of the 
        %  state equation x(k+1) = A(k)x(k) + B(k)u(k) + w(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object 
        %  K: step number
        %
        %  - Output variable(s) -
        %  VALUE: Jacobian with respect to x of the state equation. This is
        %  identical to the matrix A (at step K).
        %  
        %  - Construction -          
        %  VALUE = EVAL_FTOTJACX(OBJ,~,K,~,~) calculates the Jacobian with
        %  respect to x of the state equation. EVAL_FTOTJACX(OBJ,~,K) is
        %  also allowed.
        
            narginchk(3, 5);
            
            value=eval_A(obj,k);
            	
        end     
        
        function value=eval_htotJacX(obj,~,k,~,~)
        %EVAL_HTOTJACX Calculates Jacobian with respect to x of the 
        %  measurement equation y(k) = C(k)x(k) + D(k)u(k) + v(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object 
        %  K: step number
        %
        %  - Output variable(s) -
        %  VALUE: Jacobian with respect to x of the measurement equation. 
        %  This is identical to the matrix C (at step K).
        %  
        %  - Construction -          
        %  VALUE = EVAL_HTOTJACX(OBJ,~,K,~,~) calculates the Jacobian with
        %  respect to x of the measurement equation. EVAL_FTOTJACX(OBJ,~,K) 
        %  is also allowed.
        
            narginchk(3, 5);
            
            value=eval_C(obj,k);    
            	
        end            

        function value=eval_ftotJacW(obj,~,~,~,~) %#ok<MANU>
        %EVAL_FTOTJACW Calculates Jacobian with respect to w of the 
        %  state equation x(k+1) = A(k)x(k) + B(k)u(k) + w(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object 
        %
        %  - Output variable(s) -
        %  VALUE: Jacobian with respect to w of the state equation.
        %  
        %  - Construction -          
        %  VALUE = EVAL_FTOTJACW(OBJ,~,~,~,~) calculates the Jacobian with
        %  respect to w of the state equation. EVAL_FTOTJACW(OBJ) is also
        %  allowed.
        %
        % NOTE: EVAL_FTOTJACW always yields the identity matrix. Still, the
        % scalar 1 is always returned as value. Multiplying with the scalar 
        % 1 gives the same result as multiplying with identity matrix.
        
            narginchk(1, 5);
            
            value=1; %or, if needed change to: value=eye(obj.xSize);
            	
        end     
        
        function value=eval_htotJacV(obj,~,~,~,~) %#ok<MANU>
        %EVAL_HTOTJACV Calculates Jacobian with respect to w of the 
        %  measurement equation y(k) = C(k)x(k) + D(k)u(k) + v(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DL object 
        %
        %  - Output variable(s) -
        %  VALUE: Jacobian with respect to v of the measurement equation.
        %  
        %  - Construction -          
        %  VALUE = EVAL_HTOTJACV(OBJ,~,~,~,~) calculates the Jacobian with
        %  respect to v of the measurement equation. EVAL_HTOTJACV(OBJ) 
        %  is also allowed.
        %
        % NOTE: EVAL_HTOTJACV always yields the identity matrix. Still, the
        % scalar 1 is always returned as value. Multiplying with the scalar 
        % 1 gives the same result as multiplying with identity matrix.
        
            narginchk(1, 5);
            
            value=1; %or, if needed change to: value=eye(obj.ySize);
            	
        end            
        
        
    end %methods
    
    methods (Access = private)    
        
        function value=checkConsist(obj)
        %CHECKCONSIST Checks if object is consistent.
        
            value=0;
            if ~isempty(obj)

                if  size(obj.A,2)~=size(obj.C,2)
                	error('DA:SystemModels:ss_DL:ss_DL:colFHMismatch','A and C must have same number of columns.')
                end                
                
                if  size(obj.B,2)~=size(obj.D,2)
                	error('DA:SystemModels:ss_DL:ss_DL:colGIMismatch','B and D must have same number of columns.')
                end                                
                
                if  size(obj.A,1)~=size(obj.B,1)
                	error('DA:SystemModels:ss_DL:ss_DL:rowFGMismatch','A and B must have same number of rows.')
                end                                                
                
                if  size(obj.C,1)~=size(obj.D,1)
                	error('DA:SystemModels:ss_DL:ss_DL:rowHIMismatch','C and D must have same number of rows.')
                end          

                if obj.x0.p ~= obj.xSize
                	error('DA:SystemModels:ss_DL:ss_DL:x0StateMismatch','Size of x0 does not match with A, B and C.')
                end                    
                
                if obj.w.p ~= obj.xSize
                	error('DA:SystemModels:ss_DL:ss_DL:wInpMismatch','Size of w does not match with A, B and C.')
                end                                   
                
                if obj.v.p ~= obj.ySize
                	error('DA:SystemModels:ss_DL:ss_DL:vInpMismatch','Size of v does not match with B and D.')
                end                                                   
                
                if (length(size(obj.A))==3) && (size(obj.A,3)~=length(obj.kIndex))
                    error('DA:SystemModels:ss_DL:ss_DL:tIndexFDimMismatch','The dimensions of A and kIndex must agree')
                end             
                
                if (length(size(obj.B))==3) && (size(obj.B,3)~=length(obj.kIndex))
                    error('DA:SystemModels:ss_DL:ss_DL:tIndexGDimMismatch','The dimensions of B and kIndex must agree')
                end                       
                
                if (length(size(obj.C))==3) && (size(obj.C,3)~=length(obj.kIndex))
                    error('DA:SystemModels:ss_DL:ss_DL:tIndexHDimMismatch','The dimensions of C and kIndex must agree')
                end             
                
                if (length(size(obj.D))==3) && (size(obj.D,3)~=length(obj.kIndex))
                    error('DA:SystemModels:ss_DL:ss_DL:tIndexIDimMismatch','The dimensions of D and kIndex must agree')
                end   
                value=1;
            end
        
        end
    
    end %methods
end%classdef