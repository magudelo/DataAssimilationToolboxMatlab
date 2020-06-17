classdef ss_DNL_AN < ss_D
    %SS_DNL_AN Constructs a Discrete Non-Linear State-Space object with
    %additive noise
    %
    %      x(k+1) = f(x(k),u(k),k) + w(k)
    %      y(k)   = h(x(k),u(k),k) + v(k)    
    %
    %  SS_DNL_AN Properties:
    %     f - function handle for f(). 
    %         Syntax of function file: value = funf(x,k,u,~,Ts)
    %
    %     h - function handle for h().
    %         Syntax of function file: value = funh(x,k,u,~,Ts)
    %   
    %     fJacX - function handle for Jacobian of f() with respect to x.
    %             Syntax of function file: value = funfjacx(x,k,u,~,Ts)
    %             If not provided or [], the Jacobian is estimated numerically.
    %
    %     hJacX - function handle for Jacobian of h() with respect to x.
    %             Syntax of function file: value = funhjacx(x,k,u,~,Ts)
    %             If not provided or [], the Jacobian is estimated numerically.    
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
    %  SS_DNL_AN Construction:
    %     SS_DNL_ANMODEL=SS_DNL_AN(F,H,X0,W,V,K0,TS,TIMEUNIT,FJACX,HJACX)
    %     creates a Discrete Non-Linear State-Space object with additive
    %     noise. When an empty object is inserted for K0, TS or TIMEUNIT 
    %     their default value is used. When FJACX or HJACX are not provided
    %     their value is estimated numerically.
    %
    %     SS_DNL_ANMODEL=SS_DNL_AN(F,H,X0,W,V,K0,TS,TIMEUNIT,FJACX) estimates
    %     HJACX numerically.
    %
    %     SS_DNL_ANMODEL=SS_DNL_AN(F,H,X0,W,V,K0,TS,TIMEUNIT) estimates FJACX
    %     and HJACX numerically.
    %         
    %     SS_DNL_ANMODEL=SS_DNL_AN(F,H,X0,W,V,K0,TS) estimates FJACX and
    %     HJACX numerically and uses TIMEUNIT='seconds'.
    %      
    %     SS_DNL_ANMODEL=SS_DNL_AN(F,H,X0,W,V,K0) uses TS=1, estimates
    %     FJACX and HJACX numerically and uses TIMEUNIT='seconds'.
    %
    %     SS_DNL_ANMODEL=SS_DNL_AN(F,H,X0,W,V) uses K0=0, TS=1, estimates
    %     FJACX and HJACX numerically and uses TIMEUNIT='seconds'.
    %
    %  SS_DNL_AN Methods:
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
    %     checkConsist	- Checks if ss_DNL_AN object is properly configured.
    %
    %  See also SS_DL, SS_DNL
    
    properties (Access = public)
       
        f;          % function handle for f(). 
        fJacX;      % function handle for Jacobian of f() with respect to x.
        h;          % function handle for h(). 
        hJacX;      % function handle for Jacobian of h() with respect to x.
        x0;         % Initial state - noise model.
        w;          % Process noise - noise model.
        v;          % Measurement noise - noise model.
        k0;         % initial step number. 
        Ts;         % sample time.
        TimeUnit;   % String representing the unit of the time variable. 
           
    end
    
    methods
        
        function obj = ss_DNL_AN(f,h,x0,w,v,k0,Ts,TimeUnit,fJacX,hJacX)
        % ss_DNL_AN Constructor     
        
            narginchk(5, 10);
            ni = nargin;

            if ni >=5
                obj.f=f;
                obj.h=h;
                obj.x0=x0;
                obj.w=w;
             	obj.v=v;
            end

            if ni<6;k0=0;end;              
            if ni<7;Ts=1;end;            
            if ni<8;TimeUnit='seconds';end;        
      
            if isempty(k0);k0=0;end;               
            if isempty(Ts);Ts=1;end;            
            if isempty(TimeUnit);TimeUnit='seconds';end;   
            
            obj.k0=k0;            
            obj.Ts=Ts;
            obj.TimeUnit=TimeUnit;
            
            if ni>=9 && ~isempty(fJacX);obj.fJacX=fJacX;end;
            if ni>=10 && ~isempty(hJacX);obj.hJacX=hJacX;end;            
            
        end
        
        function obj = set.f(obj,value)
            
            value=checkArgs(value,'f');             %check f  
            obj.f=value;
            checkConsist(obj);
            
        end         
        
        function obj = set.h(obj,value)
            
            value=checkArgs(value,'h');         	%check h             
            obj.h=value;
            checkConsist(obj);
            
        end         
        
        function obj = set.fJacX(obj,value)
            
            value=checkArgs(value,'fJacX');         %check fJacX  
            obj.fJacX=value;
            checkConsist(obj);
            
        end         
        
        function obj = set.hJacX(obj,value)
            
            value=checkArgs(value,'hJacX');         %check hJacX             
            obj.hJacX=value;
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
        
        function display(obj)
        %DISPLAY Displays the object summary.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL_AN object     
        %
        %  - Construction -
        %  DISPLAY(OBJ) displays the ss_DNL_AN object summary.
        
        	disp(' ');
            disp([inputname(1),' = '])
            disp(' ');          
            if ~isempty(obj)              
                disp('Discrete Non-linear Affine Noise State-Space Model.')
                disp('---------------------------------------------------')
                disp('Properties:');
                fprintf('f (.f) \t = \t[%s] \n',f2str(obj.f));
                fprintf('h (.h) \t = \t[%s] \n',f2str(obj.h));
                fprintf('Jacobian df/dx (.fJacX) = [%s] \n',f2str(obj.fJacX));
                fprintf('Jacobian dh/dx (.hJacX) = [%s] \n',f2str(obj.hJacX));              
                fprintf('Initial State (.x0)\t = \t[%s] \n',class(obj.x0));
                fprintf('Process noise (.w) \t = \t[%s] \n',class(obj.w));
                fprintf('Measur. noise (.v) \t = \t[%s] \n',class(obj.v));
                fprintf('Initial step nr (.k0)= \t%i \n',obj.k0);                      
                fprintf('Sample Time(.Ts) \t = \t%i \n',obj.Ts);
                fprintf('Time Unit (.TimeUnit)= \t%s \n',obj.TimeUnit);
                disp(' ');
            else
                disp('Empty Discrete Non-linear Affine Noise State-Space Model.')
            end                
        end        
        
        function value=isempty(obj)
        %ISEMPTY Check if object is empty.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL_AN object   
        %
        %  - Output variable(s) -
        %  VALUE: indicates if object is empty (=1) or not empty (=0).
        %  
        %  - Construction -          
        %  VALUE = ISEMPTY(OBJ) checks if object is empty or not. The
        %  object is considered empty when f, h, x0, w, v, Ts or TimeUnit
        %  are empty.
        
            if (~isempty(obj.f)&& ~isempty(obj.h) && ~isempty(obj.x0)&& ~isempty(obj.w)&&...
                ~isempty(obj.v)&& ~isempty(obj.Ts)&& ~isempty(obj.TimeUnit))
                value=0;
            else
                value=1;
            end

        end                
        
        function value=eval_ftot(obj,x,k,u,w)
        %EVAL_FTOT Calculates x(k+1) of state equation.
        %   x(k+1) = f(x(k),u(k),k) + w(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL_AN object 
        %  X: state at step K
        %  K: step number
        %  U: input at step K
        %  W: process noise at step K. (W=[] takes sample from distribution, 
        %     W=0 calculates f(x(k),u(k),k) )
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
            
            value=obj.f(x,k,u,[],obj.Ts);   
            
            if isempty(w)           	%if w=[]: take sample  
            	w=sample(obj.w,1,k);
                value = value + w;
            elseif ~isequal(w,0)        %if w=0: add nothing
                value = value + w; 
            end
            	
        end          
        
        function value=eval_htot(obj,x,k,u,v)
        %EVAL_HTOT Calculates x(k+1) of measurement equation.
        %   y(k) = h(x(k),u(k),k) + v(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL_AN object 
        %  X: state at step K
        %  K: step number
        %  U: input at step K
        %  V: measurement noise at step K. (V=[] takes sample from  
        %     distribution, V=0 calculates h(x(k),u(k),k) )
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
            
            value=obj.h(x,k,u,[],obj.Ts);  
            
            if isempty(v)           	%if v=[]: take sample  
                v=sample(obj.v,1,k);
                value = value + v;
            elseif ~isequal(v,0)        %if v=0: add nothing
                value = value + v; 
            end
            	
        end           

        function value=eval_ftotJacX(obj,x,k,u,~)
        %EVAL_FTOTJACX Calculates Jacobian with respect to x of the 
        %  state equation x(k+1) = f(x(k),u(k),k) + w(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL_AN object 
        %  X: state at step K
        %  K: step number        
        %  U: input at step K
        %
        %  - Output variable(s) -
        %  VALUE: Jacobian with respect to x of the state equation. When
        %  no function handle for the Jacobian is provided, the Jacobian is
        %  numerically estimated.
        %  
        %  - Construction -          
        %  VALUE = EVAL_FTOTJACX(OBJ,X,K,U,~) calculates the Jacobian with
        %  respect to x of the state equation. EVAL_FTOTJACX(OBJ,X,K,U) is
        %  also allowed.
        
            narginchk(4, 5);
            
            if ~isempty(obj.fJacX)
                value=obj.fJacX(x,k,u,[],obj.Ts); 
            else
                value = jacest(obj.f,'x',x,k,u,[],obj.Ts);
            end
        end     
        
        function value=eval_htotJacX(obj,x,k,u,~)
        %EVAL_HTOTJACX Calculates Jacobian with respect to x of the 
        %  measurement equation y(k) = h(x(k),u(k),k) + v(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL_AN object 
        %  X: state at step K
        %  K: step number        
        %  U: input at step K
        %
        %  - Output variable(s) -
        %  VALUE: Jacobian with respect to x of the measurement equation. 
        %  When no function handle for the Jacobian is provided, the 
        %  Jacobian is numerically estimated.
        %  
        %  - Construction -          
        %  VALUE = EVAL_HTOTJACX(OBJ,X,K,U,~) calculates the Jacobian with
        %  respect to x of the measurement equation. EVAL_HTOTJACX(OBJ,X,K,U)
        %  is also allowed.
        
            narginchk(4, 5);
            
            if ~isempty(obj.hJacX)
                value=obj.hJacX(x,k,u,[],obj.Ts); 
            else
                value = jacest(obj.h,'x',x,k,u,[],obj.Ts);
            end
            	
        end            
        
        function value=eval_ftotJacW(obj,~,~,~,~) %#ok<MANU>
        %EVAL_FTOTJACW Calculates Jacobian with respect to w of the 
        %  state equation x(k+1) = f(x(k),u(k),k) + w(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL_AN object 
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
            
            value=1;
            	
        end     
        
        function value=eval_htotJacV(obj,~,~,~,~) %#ok<MANU>
        %EVAL_HTOTJACV Calculates Jacobian with respect to w of the 
        %  measurement equation y(k) = h(x(k),u(k),k) + v(k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL_AN object 
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
            
            value=1;
            	
        end            
        
        function value=checkConsist(obj,varargin)
        %CHECKCONSIST Checks if ss_DNL_AN object is properly configured.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL_AN object 
        %  USIZE: the input size
        %  YSIZE: the output size
        %
        %  - Output variable(s) -
        %  VALUE: indicates if object is configured properly. 
        %  (1=consistent, 0= not consistent)
        %  
        %  - Construction -          
        %  VALUE = CHECKCONSIST(OBJ,USIZE,YSIZE) checks if object is 
        %  configured properly.
        
            ni=nargin;
            value=0;
            if ni>1
                xSize = obj.x0.p;
                uSize = varargin{1};
                ySize = varargin{2};                
                wSize = obj.w.p;
                vSize = obj.v.p;
                
                if ~isempty(obj)

                    if wSize~=xSize
                        error('DA:SystemModels:ss_DNL_AN:ss_DNL_AN:wErr','Size of w does not match with amount of system states.')
                    end
                    
                    if vSize~=ySize
                        error('DA:SystemModels:ss_DNL_AN:ss_DNL_AN:vErr','Size of v does not match with amount of system outputs.')
                    end                    
                    
                    try
                        val=obj.f(zeros(xSize,1),obj.k0,zeros(uSize,1),[],obj.Ts);
                        if length(val)~=xSize
                            error('DA:SystemModels:ss_DNL_AN:ss_DNL_AN:fErr','Error in function f: function output does not match with amount of system states.')
                        end                     
                    catch err
                        error('DA:SystemModels:ss_DNL_AN:ss_DNL_AN:fErr','Error in function f: %s',err.message)
                    end

                    try                
                        val=obj.h(zeros(xSize,1),obj.k0,zeros(uSize,1),[],obj.Ts);
                        if length(val)~=ySize
                            error('DA:SystemModels:ss_DNL_AN:ss_DNL_AN:hErr','Error in function h: function output does not match with amount of system outputs.')
                        end                          
                    catch err
                        error('DA:SystemModels:ss_DNL_AN:ss_DNL_AN:hErr','Error in function h: %s',err.message)
                    end    

                    if ~isempty(obj.fJacX)
                        try
                            val=obj.fJacX(zeros(xSize,1),obj.k0,zeros(uSize,1),[],obj.Ts);
                            if size(val,1)~= xSize || size(val,2)~= xSize
                                error('DA:SystemModels:ss_DNL_AN:ss_DNL_AN:fJacXErr','Error in function fJacX: Jacobian dimensions do not match with amount of states.')
                            end                                 
                        catch err
                            error('DA:SystemModels:ss_DNL_AN:ss_DNL_AN:fJacXErr','Error in function fJacX: %s',err.message)
                        end
                    end

                    if ~isempty(obj.hJacX)
                        try
                            val=obj.hJacX(zeros(xSize,1),obj.k0,zeros(uSize,1),[],obj.Ts);
                            if size(val,1)~= ySize || size(val,2)~= xSize
                                error('DA:SystemModels:ss_DNL_AN:ss_DNL_AN:hJacXErr','Error in function hJacX: Jacobian dimensions do not match with amount of states and outputs.')
                            end                              
                        catch err
                            error('DA:SystemModels:ss_DNL_AN:ss_DNL_AN:hJacXErr','Error in function hJacX: %s',err.message)
                        end
                    end       
                    
                    value=1;
                end
            end
        end        
        
    end %methods
    
    methods (Access = private)    
    
    end %methods
end%classdef