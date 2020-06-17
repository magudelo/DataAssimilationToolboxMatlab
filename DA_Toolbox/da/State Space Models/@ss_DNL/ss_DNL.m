classdef ss_DNL < ss_D
    %SS_DNL Constructs a Discrete Non-Linear State-Space object 
    %
    %      x(k+1) = f(x(k),u(k),w(k),k) 
    %      y(k)   = h(x(k),u(k),v(k),k)    
    %
    %  SS_DNL Properties:
    %     f - function handle for f(). 
    %         Syntax of function file: value = funf(x,k,u,w,Ts)
    %
    %     h - function handle for h().
    %         Syntax of function file: value = funh(x,k,u,v,Ts)
    %   
    %     fJacX - function handle for Jacobian of f() with respect to x.
    %             Syntax of function file: value = funfjacx(x,k,u,w,Ts)
    %             If not provided or [], the Jacobian is estimated numerically.
    %
    %     fJacW - function handle for Jacobian of f() with respect to w.
    %             Syntax of function file: value = funfjacw(x,k,u,w,Ts)
    %             If not provided or [], the Jacobian is estimated numerically.
    %    
    %     hJacX - function handle for Jacobian of h() with respect to x.
    %             Syntax of function file: value = funhjacx(x,k,u,v,Ts)
    %             If not provided or [], the Jacobian is estimated numerically.    
    %   
    %     hJacV - function handle for Jacobian of h() with respect to v.
    %             Syntax of function file: value = funhjacx(x,k,u,v,Ts)
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
    %  SS_DNL Construction:
    %     SS_DNLMODEL=SS_DNL(F,H,X0,W,V,K0,TS,TIMEUNIT,FJACX,FJACW,HJACX,HJACV)
    %     creates a Discrete Non-Linear State-Space object. When an empty
    %     object is inserted for K0, TS or TIMEUNIT their default value
    %     is used. When FJACX, FJACW, HJACX or HJACV are not provided
    %     their value is estimated numerically.
    %
    %     SS_DNLMODEL=SS_DNL(F,H,X0,W,V,K0,TS,TIMEUNIT,FJACX,FJACW,HJACX) 
    %     estimates HJACV numerically.
    %
    %     SS_DNLMODEL=SS_DNL(F,H,X0,W,V,K0,TS,TIMEUNIT,FJACX,FJACW) 
    %     estimates HJACX, HJACV numerically.
    %        
    %     SS_DNLMODEL=SS_DNL(F,H,X0,W,V,K0,TS,TIMEUNIT,FJACX) 
    %     estimates HJACX, HJACV and FJACW numerically.
    %     
    %     SS_DNLMODEL=SS_DNL(F,H,X0,W,V,K0,TS,TIMEUNIT) 
    %     estimates HJACX, HJACV, FJACX and FJACW numerically.
    %        
    %     SS_DNLMODEL=SS_DNL(F,G,H,I,X0,W,V,K0,TS) estimates the Jacobians
    %     numerically and uses TIMEUNIT='seconds'.
    %      
    %     SS_DNLMODEL=SS_DNL(F,G,H,I,X0,W,V,K0) uses TS=1,
    %     TIMEUNIT='seconds' and estimates the Jacobians numerically.
    %
    %     SS_DNLMODEL=SS_DNL(F,G,H,I,X0,W,V) uses K0=0, TS=1,
    %     TIMEUNIT='seconds' and estimates the Jacobians numerically.
    %
    %  SS_DNL Methods:
    %     eval_ftot - Calculates x(k+1)=f() of state equation.
    %     eval_htot - Calculates y(k)=h() of measurement equation.	
    %     eval_ftotJacX - Calculates Jacobian with respect to x of the 
    %                     state equation f().
    %     eval_htotJacX - Calculates Jacobian with respect to x of the 
    %                     measurement equation h().  
    %     eval_ftotJacW - Calculates Jacobian with respect to w of the 
    %                     state equation f().
    %     eval_htotJacV - Calculates Jacobian with respect to v of the 
    %                     measurement equation h().
    %     checkConsist	- Checks if SS_DNL object is properly configured.
    %
    %  See also SS_DL, SS_DNL_AN
    
    properties (Access = public)
       
        f;       	% function handle for f(). 
        fJacX;      % function handle for Jacobian of f() with respect to x.
        fJacW;      % function handle for Jacobian of f() with respect to w.
        h;          % function handle for h(). 
        hJacX;      % function handle for Jacobian of h() with respect to x.
        hJacV;      % function handle for Jacobian of h() with respect to v.
        x0;         % Initial state - noise model.
        w;          % Process noise - noise model.
        v;          % Measurement noise - noise model.
        k0;         % initial step number. 
        Ts;         % sample time.
        TimeUnit;   % String representing the unit of the time variable. 
           
    end
    
    methods
        
        function obj = ss_DNL(f,h,x0,w,v,k0,Ts,TimeUnit,fJacX,fJacW,hJacX,hJacV)
        % ss_DNL Constructor  
        
            narginchk(5, 12);
            ni = nargin;

            if ni>=5
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
            
            if ni>=8 && ~isempty(fJacX);obj.fJacX=fJacX;end;
            if ni>=9 && ~isempty(fJacW);obj.fJacW=fJacW;end;            
            if ni>=10 && ~isempty(hJacX);obj.hJacX=hJacX;end;
            if ni>=11 && ~isempty(hJacV);obj.hJacV=hJacV;end; 
            
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
        
        function obj = set.fJacW(obj,value)
            
            value=checkArgs(value,'fJacW');         %check fJacX  
            obj.fJacW=value;
            checkConsist(obj);
            
        end             
        
        function obj = set.hJacX(obj,value)
            
            value=checkArgs(value,'hJacX');         %check hJacX             
            obj.hJacX=value;
            checkConsist(obj);
            
        end          
        
        function obj = set.hJacV(obj,value)
            
            value=checkArgs(value,'hJacV');         %check hJacX             
            obj.hJacV=value;
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
        %  OBJ: a ss_DNL object     
        %
        %  - Construction -
        %  DISPLAY(OBJ) displays the ss_DNL object summary.
        
        	disp(' ');
            disp([inputname(1),' = '])
            disp(' ');          
            if ~isempty(obj)              
                disp('Discrete Non-Linear State-Space Model.')
                disp('--------------------------------------')
                disp('Properties:');
                fprintf('f (.f) \t = \t[%s] \n',f2str(obj.f));
                fprintf('h (.h) \t = \t[%s] \n',f2str(obj.h));
                fprintf('Jacobian df/dx (.fJacX) = [%s] \n',f2str(obj.fJacX));
                fprintf('Jacobian df/dw (.fJacW) = [%s] \n',f2str(obj.fJacW));                
                fprintf('Jacobian dh/dx (.hJacX) = [%s] \n',f2str(obj.hJacX));              
                fprintf('Jacobian dh/dv (.hJacV) = [%s] \n',f2str(obj.hJacV));                
                fprintf('Initial State (.x0)\t = \t[%s] \n',class(obj.x0));
                fprintf('Process noise (.w) \t = \t[%s] \n',class(obj.w));
                fprintf('Measur. noise (.v) \t = \t[%s] \n',class(obj.v));
                fprintf('Initial step nr (.k0)= \t%i \n',obj.k0);                     
                fprintf('Sample Time(.Ts) \t = \t%i \n',obj.Ts);
                fprintf('Time Unit (.TimeUnit)= \t%s \n',obj.TimeUnit);
                disp(' ');
            else
                disp('Empty Discrete Non-Linear State-Space Model.')
            end                
        end        
        
        function value=isempty(obj)
        %ISEMPTY Check if object is empty.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL object   
        %
        %  - Output variable(s) -
        %  VALUE: indicates if object is empty (=1) or not empty (=0).
        %  
        %  - Construction -          
        %  VALUE = ISEMPTY(OBJ) checks if object is empty or not. The
        %  object is considered empty when f, h, x0, w, v, Ts or TimeUnit
        %  are empty.
        
            if (~isempty(obj.f)&& ~isempty(obj.h)&& ~isempty(obj.x0)&& ~isempty(obj.w)&&...
                ~isempty(obj.v)&& ~isempty(obj.Ts)&& ~isempty(obj.TimeUnit))
                value=0;
            else
                value=1;
            end

        end                
        
        function value=eval_ftot(obj,x,k,u,w)
        %EVAL_FTOT Calculates x(k+1) of state equation.
        %   x(k+1) = f(x(k),u(k),w(k),k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL object 
        %  X: state at step K
        %  K: step number
        %  U: input at step K
        %  W: process noise at step K. (W=[] takes sample from distribution
        %  , W=0 can be entered for the zero vector)
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

            if isempty(w)                   %if w=[]: take sample           
                w=sample(obj.w,1,k);
            elseif isequal(w,0)             %if w=0: add nothing
                w =zeros(obj.w.p,1);                
            end
            
            value=obj.f(x,k,u,w,obj.Ts);   
            	
        end          
        
        function value=eval_htot(obj,x,k,u,v)
        %EVAL_HTOT Calculates x(k+1) of measurement equation.
        %   y(k) = h(x(k),u(k),v(k),k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL object 
        %  X: state at step K
        %  K: step number
        %  U: input at step K
        %  V: measurement noise at step K. (V=[] takes sample from  
        %     distribution, V=0 can be entered for the zero vector)
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
            
            if isempty(v)               %if v=[]: take sample
                v=sample(obj.v,1,k);
            elseif isequal(v,0)         %if v=0: add nothing
                v =zeros(obj.v.p,1);                
            end
            
            value=obj.h(x,k,u,v,obj.Ts);  
            	
        end           

         function value=eval_ftotJacX(obj,x,k,u,w)
        %EVAL_FTOTJACX Calculates Jacobian with respect to x of the 
        %  state equation x(k+1) = f(x(k),u(k),w(k),k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL object 
        %  X: state at step K
        %  K: step number        
        %  U: input at step K
        %  W: process noise at step K. (W=[] takes sample from distribution
        %  , W=0 can be entered for the zero vector)        
        %
        %  - Output variable(s) -
        %  VALUE: Jacobian with respect to x of the state equation. When
        %  no function handle for the Jacobian is provided, the Jacobian is
        %  numerically estimated.
        %  
        %  - Construction -          
        %  VALUE = EVAL_FTOTJACX(OBJ,X,K,U,W) calculates the Jacobian with
        %  respect to x of the state equation. 
        %
        %  VALUE = EVAL_FTOTJACX(OBJ,X,K,U) calculates the Jacobian with
        %  respect to x of the state equation with process noise W a 
        %  sample from its distribution.       
        
            narginchk(4, 5);
            ni = nargin;  
            
            if ni~=5;w=[];end;

            if isempty(w)               %if w=[]: take sample
                w=sample(obj.w,1,k);
            elseif isequal(w,0)         %if w=0: add nothing
                w =zeros(obj.w.p,1);                
            end            
            
            if  ~isempty(obj.fJacX)
                value=obj.fJacX(x,k,u,w,obj.Ts); 
            else
                value = jacest(obj.f,'x',x,k,u,w,obj.Ts);
            end
        end     
        
        function value=eval_htotJacX(obj,x,k,u,v)
        %EVAL_HTOTJACX Calculates Jacobian with respect to x of the 
        %  measurement equation y(k) = h(x(k),u(k),v(k),k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL object 
        %  X: state at step K
        %  K: step number        
        %  U: input at step K
        %  V: measurement noise at step K. (V=[] takes sample from  
        %     distribution, V=0 can be entered for the zero vector)
        %
        %  - Output variable(s) -
        %  VALUE: Jacobian with respect to x of the measurement equation. 
        %  When no function handle for the Jacobian is provided, the 
        %  Jacobian is numerically estimated.
        %  
        %  - Construction -          
        %  VALUE = EVAL_HTOTJACX(OBJ,X,K,U,V) calculates the Jacobian with
        %  respect to x of the measurement equation. 
        %
        %  VALUE = EVAL_HTOTJACX(OBJ,X,K,U) calculates the Jacobian with
        %  respect to x of the measurement equation with measurement noise 
        %  V a sample from its distribution.         
        
            narginchk(4, 5);
            ni = nargin;  
            
            if ni~=5;v=[];end;

            if isempty(v)               %if v=[]: take sample
                v=sample(obj.w,1,k);
            elseif isequal(v,0)         %if v=0: add nothing
                v =zeros(obj.w.p,1);                
            end                  
            
            if ~isequal(obj.hJacX,0) && ~isempty(obj.hJacX)
                value=obj.hJacX(x,k,u,v,obj.Ts); 
            else
                value = jacest(obj.h,'x',x,k,u,v,obj.Ts);
            end
            	
        end         

         function value=eval_ftotJacW(obj,x,k,u,w)
        %EVAL_FTOTJACW Calculates Jacobian with respect to w of the 
        %  state equation x(k+1) = f(x(k),u(k),w(k),k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL object 
        %  X: state at step K
        %  K: step number        
        %  U: input at step K
        %  W: process noise at step K. (W=[] takes sample from distribution
        %  , W=0 can be entered for the zero vector)        
        %
        %  - Output variable(s) -
        %  VALUE: Jacobian with respect to w of the state equation. When
        %  no function handle for the Jacobian is provided, the Jacobian is
        %  numerically estimated.
        %  
        %  - Construction -          
        %  VALUE = EVAL_FTOTJACW(OBJ,X,K,U,W) calculates the Jacobian with
        %  respect to w of the state equation. 
        %
        %  VALUE = EVAL_FTOTJACW(OBJ,X,K,U) calculates the Jacobian with
        %  respect to w of the state equation with process noise W a 
        %  sample from its distribution.  
        
            narginchk(4, 5);
            ni = nargin;  
            
            if ni~=5;w=[];end;

            if isempty(w)               %if w=[]: take sample
                w=sample(obj.w,1,k);
            elseif isequal(w,0)         %if w=0: add nothing
                w =zeros(obj.w.p,1);                
            end            
            
            if ~isequal(obj.fJacW,0) && ~isempty(obj.fJacW)
                value=obj.fJacW(x,k,u,w,obj.Ts); 
            else
                value = jacest(obj.f,'w',x,k,u,w,obj.Ts);
            end
        end     
        
        function value=eval_htotJacV(obj,x,k,u,v)
        %EVAL_HTOTJACV Calculates Jacobian with respect to v of the 
        %  measurement equation y(k) = h(x(k),u(k),v(k),k)
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL object 
        %  X: state at step K
        %  K: step number        
        %  U: input at step K
        %  V: measurement noise at step K. (V=[] takes sample from  
        %     distribution, V=0 can be entered for the zero vector)
        %
        %  - Output variable(s) -
        %  VALUE: Jacobian with respect to v of the measurement equation. 
        %  When no function handle for the Jacobian is provided, the 
        %  Jacobian is numerically estimated.
        %  
        %  - Construction -          
        %  VALUE = EVAL_HTOTJACV(OBJ,X,K,U,V) calculates the Jacobian with
        %  respect to v of the measurement equation. 
        %
        %  VALUE = EVAL_HTOTJACV(OBJ,X,K,U) calculates the Jacobian with
        %  respect to v of the measurement equation with measurement noise 
        %  V a sample from its distribution.    
        
            narginchk(4, 5);
            ni = nargin;  
            
            if ni~=5;v=[];end;

            if isempty(v)               %if v=[]: take sample
                v=sample(obj.v,1,k);
            elseif isequal(v,0)         %if v=0: add nothing
                v =zeros(obj.v.p,1);                
            end                  
            
            if ~isequal(obj.hJacV,0) && ~isempty(obj.hJacV)
                value=obj.hJacV(x,k,u,v,obj.Ts); 
            else
                value = jacest(obj.h,'v',x,k,u,v,obj.Ts);
            end
            	
        end             
        
        function value=checkConsist(obj,varargin)
        %CHECKCONSIST Checks if ss_DNL object is properly configured.
        %
        %  - Input variable(s) -
        %  OBJ: a ss_DNL object 
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

                    try
                        val=obj.f(zeros(xSize,1),obj.k0,zeros(uSize,1),zeros(wSize,1),obj.Ts);
                        if length(val)~=xSize
                            error('DA:SystemModels:ss_DNL:ss_DNL:fErr','Error in function f: function output does not match with amount of system states.')
                        end                          
                    catch err
                        error('DA:SystemModels:ss_DNL:ss_DNL:fErr','Error in function f: %s',err.message)
                    end

                    try                
                        val=obj.h(zeros(xSize,1),obj.k0,zeros(uSize,1),zeros(vSize,1),obj.Ts);
                        if length(val)~=ySize
                            error('DA:SystemModels:ss_DNL:ss_DNL:hErr','Error in function h: function output does not match with amount of system outputs.')
                        end                            
                    catch err
                     	error('DA:SystemModels:ss_DNL:ss_DNL:hErr','Error in function h: %s',err.message)
                    end    

                    if ~isempty(obj.fJacX)
                        try
                            val=obj.fJacX(zeros(xSize,1),obj.k0,zeros(uSize,1),obj.Ts);
                            if size(val,1)~= xSize || size(val,2)~= xSize
                                error('DA:SystemModels:ss_DNL:ss_DNL:fJacXErr','Error in function fJacX: Jacobian dimensions do not match with amount of states.')
                            end                                
                        catch err
                            error('DA:SystemModels:ss_DNL:ss_DNL:fJacXErr','Error in function fJacX: %s',err.message)
                        end
                    end

                    if ~isempty(obj.hJacX)
                        try
                            val=obj.hJacX(zeros(xSize,1),obj.k0,zeros(uSize,1),obj.Ts);
                            if size(val,1)~= ySize || size(val,2)~= xSize
                                error('DA:SystemModels:ss_DNL:ss_DNL:hJacXErr','Error in function hJacX: Jacobian dimensions do not match with amount of states and outputs.')
                            end                                  
                        catch err
                            error('DA:SystemModels:ss_DNL:ss_DNL:hJacXErr','Error in function hJacX: %s',err.message)
                        end
                    end        

                    if ~isempty(obj.fJacW)
                        try
                            val=obj.fJacW(zeros(xSize,1),obj.k0,zeros(uSize,1),obj.Ts);
                            if size(val,1)~= xSize || size(val,2)~= wSize
                                error('DA:SystemModels:ss_DNL:ss_DNL:fJacWErr','Error in function fJacW: Jacobian dimensions do not match with amount of states and size of w.')
                            end                               
                        catch err
                            error('DA:SystemModels:ss_DNL:ss_DNL:fJacWErr','Error in function fJacW: %s',err.message)
                        end
                    end

                    if ~isempty(obj.hJacV)
                        try
                            val=obj.hJacV(zeros(xSize,1),obj.k0,zeros(uSize,1),obj.Ts);
                            if size(val,1)~= ySize || size(val,2)~= vSize
                                error('DA:SystemModels:ss_DNL:ss_DNL:hJacVErr','Error in function hJacV: Jacobian dimensions do not match with amount of outputs and size of v.')
                            end                               
                        catch err
                            error('DA:SystemModels:ss_DNL:ss_DNL:hJacVErr','Error in function hJacV: %s',err.message)
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