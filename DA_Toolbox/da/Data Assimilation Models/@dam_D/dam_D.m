classdef dam_D
    %DAM_D Constructs a discrete data assimilation model
    %
    %  DAM_D Properties:
    %
    %     alg - string that defines the algorithm used to generate the
    %     data. e.g. KF, EKF, UKF,... (read-only) 
    %
    %     x - True states. Matrix of true states where each column is 
    %     a true state. (Can be empty []) (read-only) 
    %                
    %     xf - states forecast. Matrix of assimilated forecast
    %     states where each column is a forecast state. (Can be empty [])
    %     (read-only) 
    %
    %     xa - states analysis. Matrix of assimilated analysis
    %     states where each column is a analysis state. (must be non-empty) 
    %     (read-only) 
    %
    %     Pf - (co)variance forecast. 3D array of assimilated 
    %     forecast covariance matrices (covFull=1) or matrix of assimilated 
    %     forecast variances where each column is a variance (covFull=0).
    %     (Can be empty []) (read-only) 
    %
    %     Pa - (co)variance analysis. 3D array of assimilated 
    %     analysis covariance matrices (covFull=1) or matrix of assimilated 
    %     analysis variances where each column is a variance (covFull=0).    
    %     (Can be empty []) (read-only)  
    %
    %     u - Inputs. Matrix of inputs where each column is an input.
    %     (Can be empty []) (read-only) 
    %
    %     y - Measurements. Matrix of measurements where each column is
    %     a measurement. (Can be empty []) (read-only) 
    %
    %     k - step number array. Each element contains the step number
    %     corresponding to each column of xf, xa, Pf, Pa, u, y... 
    %     For example: the second element of k corresponds to the step 
    %     number of column two. (must be non-empty) 
    %
    %     Ts - sample time. (default is 1) (read-only)   
    %
    %     covFull - Variable that indicates whether the full covariance
    %     matrices are saved (=1) or the only the variances (=0). 
    %     Is required as input argument if Pf or Pa are non-empty.
    %     (access=private)
    %
    %     TimeUnit - String representing the unit of the time variable. 
    %       TimeUnit can take the following values: 'nanoseconds',
    %       'microseconds', 'milliseconds', 'seconds', 'minutes', 'hours',
    %       'days', 'weeks', 'months', 'years'. (read-only) 
    %       Default: 'seconds'     
    %
    %     l - Total amount of iterations.  (read-only) 
    %    
    %  DAM_D Construction:
    %     OBJ = DAM_D(ALG,X,XF,XA,PF,PA,U,Y,K0,COVFULL,TS,TIMEUNIT)  
    %     creates an object OBJ representing a discrete data assimilation
    %     model.
    %
    %     OBJ = DAM_D(ALG,X,XF,XA,PF,PA,U,Y,K0,COVFULL,TS) creates an 
    %     object OBJ representing a discrete data assimilation model with 
    %     TIMEUNIT='seconds'.
    %
    %     OBJ = DAM_D(ALG,X,XF,XA,PF,PA,U,Y,K0,COVFULL) creates an object
    %     OBJ representing a discrete data assimilation model with TS=1 and
    %     TIMEUNIT='seconds'.
    %    
    %     OBJ = DAM_D(ALG,X,XF,XA,[],[],U,Y,K0) creates an object OBJ 
    %     representing a discrete data assimilation model with TS=1 and
    %     TIMEUNIT='seconds'. Since PA and PF are empty COVFULL is not 
    %     required.
    %        
    %  DAM_D Methods:
    %     plot - Overloaded plot method to plot dam_D objects. 
    %     rmse - Calculates the RMSE between x and xa of dam_D objects. 
    %     mse - Calculates the MSE between x and xa of dam_D objects.    
    % 
    
	properties (SetAccess = private)

        alg;            % Algorithm used to create data.
        x;              % True states matrix.
        xf;             % Assimilated states forecast matrix
        xa;             % Assimilated states analysis matrix (must be non-empty)         
        Pf;             % Assimilated (co)variance forecast 3D array of matrices.
        Pa;             % Assimilated (co)variance analysis 3D array of matrices.     
        u;              % Inputs matrix
        y;              % Measurements 
        k;              % Step number array (must be non-empty)      
        Ts;             % Sample time
        TimeUnit;       % String representing the unit of the time variable.       
        
    end
    
    properties (Access = private)
       covFull;         %full covariance matrices are saved (=1) or the only the variances (=0).
    end

	properties (Dependent = true, SetAccess = private)
        l;              % Total amount of iterations.
    end       
    
    methods
        
        function obj = dam_D(alg,x,xf,xa,Pf,Pa,u,y,k,covFull,Ts,TimeUnit)
        % DAM_D Constructor 
        
            narginchk(9,12);
            
            ni=nargin;
            
            obj.alg=alg;
            obj.x=x;
            obj.xf=xf;      
            obj.xa=xa;  
            obj.Pf=Pf;  
            obj.Pa=Pa;                  
            obj.u=u;
            obj.y=y;
            obj.k=k;           

            if ni<10;covFull=-1;end;            
            if ni<11;Ts=1;end;      
            if ni<12;TimeUnit='seconds';end;    
            if isempty(covFull);covFull=-1;end;              
            if isempty(Ts);Ts=1;end; 
            if isempty(TimeUnit);TimeUnit='seconds';end;     
           
            obj.covFull=covFull;            
            obj.Ts=Ts;
            obj.TimeUnit=TimeUnit;

        end            
        
        function obj = set.alg(obj,value)
            
            value=checkArgs(value,'alg');     	%check alg                
            obj.alg=value;
            checkConsist(obj);            
        end              
        
        function obj = set.x(obj,value)
            
            value=checkArgs(value,'x');     	%check x                
            obj.x=value;
            checkConsist(obj);            
        end             
        
        function obj = set.xf(obj,value)
            
            value=checkArgs(value,'xf');     	%check xf                
            obj.xf=value;
            checkConsist(obj);            
        end          
        
        function obj = set.xa(obj,value)
            
            value=checkArgs(value,'xa');     	%check xa                
            obj.xa=value;
            checkConsist(obj);            
        end                  
       
        function obj = set.Pf(obj,value)
            
            value=checkArgs(value,'Pf');     	%check Pf                
            obj.Pf=value;
            checkConsist(obj);            
        end          
        
        function obj = set.Pa(obj,value)
            
            value=checkArgs(value,'Pa');     	%check Pa                
            obj.Pa=value;
            checkConsist(obj);            
        end                   
        
        function obj = set.u(obj,value)
            
            value=checkArgs(value,'u');         %check u                    
            obj.u=value;
            checkConsist(obj);            
        end
        
        function obj = set.y(obj,value)
            
            value=checkArgs(value,'y');         %check y                    
            obj.y=value;
            checkConsist(obj);            
        end       
        
        function obj = set.k(obj,value)
            
            value=checkArgs(value,'k');         %check k                   
            obj.k=value;
            checkConsist(obj);            
        end             
        
        function obj = set.Ts(obj,value)
            
            value=checkArgs(value,'Ts');        %check Ts             
            obj.Ts=value;
            checkConsist(obj);
            
        end    
        
        function obj = set.covFull(obj,value)
            
            value=checkArgs(value,'covFull'); 	%check covFull             
            obj.covFull=value;
            checkConsist(obj);
            
        end            
        
        function obj = set.TimeUnit(obj,value)
            
            value=checkArgs(value,'TimeUnit'); 	%check TimeUnit                   
            obj.TimeUnit=value;
            checkConsist(obj);
            
        end                
        
        function value = get.l(obj)
            
            n=size(obj.Pf);
            if length(n)>2
                sizePf=n(3);
            else
                sizePf=n(2);
            end
            
            n=size(obj.Pa);
            if length(n)>2
                sizePa=n(3);
            else
                sizePa=n(2);
            end            
            
            value = max([size(obj.xf,2) size(obj.xa,2) sizePf sizePa size(obj.u,2) size(obj.y,2)]);
            
        end                
        
        function display(obj)
        %DISPLAY Displays the object summary.
        %
        %  - Input variable(s) -
        %  OBJ: a dam_D object     
        %
        %  - Construction -
        %  DISPLAY(OBJ) displays the dam_D object summary.
        
        	disp(' ');
            disp([inputname(1),' = '])
            disp(' ');          
            if ~isempty(obj)              
                disp('Discrete Assimilation Data Object.')
                disp('----------------------------------')
                disp('Properties:');
                if ~isempty(obj.x);fprintf('True state (.x) \t = \t[%i x %i %s] \n',size(obj.x,1),size(obj.x,2),class(obj.x));end;                
                fprintf('State analysis (.xa) \t = \t[%i x %i %s] \n',size(obj.xa,1),size(obj.xa,2),class(obj.xa));
                if ~isempty(obj.xf);fprintf('State forecast (.xf) \t = \t[%i x %i %s] \n',size(obj.xf,1),size(obj.xf,2),class(obj.xf));end;
                if length(size(obj.Pf))==3
                    if ~isempty(obj.Pa);fprintf('Covariance analysis(.Pa) = \t[%i x %i x %i %s] \n',size(obj.Pa,1),size(obj.Pa,2),size(obj.Pa,3),class(obj.Pa));end;
                    if ~isempty(obj.Pf);fprintf('Covariance forecast(.Pf) = \t[%i x %i x %i %s] \n',size(obj.Pf,1),size(obj.Pf,2),size(obj.Pa,3),class(obj.Pf));end;                      
                else
                    if ~isempty(obj.Pa);fprintf('Covariance analysis(.Pa) = \t[%i x %i %s] \n',size(obj.Pa,1),size(obj.Pa,2),class(obj.Pa));end;
                    if ~isempty(obj.Pf);fprintf('Covariance forecast(.Pf) = \t[%i x %i %s] \n',size(obj.Pf,1),size(obj.Pf,2),class(obj.Pf));  end;                                        
                end
                if ~isempty(obj.u);fprintf('Inputs (.u) \t \t \t = \t[%i x %i %s] \n',size(obj.u,1),size(obj.u,2),class(obj.u));end;
                if ~isempty(obj.y);fprintf('Measurements (.y) \t \t = \t[%i x %i %s] \n',size(obj.y,1),size(obj.y,2),class(obj.y));end;
                fprintf('Step nr index (.k) \t = \t[%i x %i %s] \n',size(obj.k,1),size(obj.k,2),class(obj.k));
                fprintf('Samples (.l) \t \t = \t%i \n',obj.l);                 
                fprintf('Sample Time(.Ts) \t = \t%i \n',obj.Ts);                
                fprintf('Time Unit(.TimeUnit) = \t%s \n',obj.TimeUnit);                
            else
                disp('Empty Discrete Assimilation Data Object.')                
            end
        end   
        
        function value=isempty(obj)
        %ISEMPTY Check if object is empty.
        %
        %  - Input variable(s) -
        %  OBJ: a dam_D object   
        %
        %  - Output variable(s) -
        %  VALUE: indicates if object is empty (=1) or not empty (=0).
        %  
        %  - Construction -          
        %  VALUE = ISEMPTY(OBJ) checks if object is empty or not. The
        %  object is considered empty when xa or k are empty.
        %  
        
            if (~isempty(obj.xa) && ~isempty(obj.k) && ~isempty(obj.covFull))
                value=0;
            else
                value=1;
            end

        end          
        
        function rmseval=rmse(obj,nrs)
        %RMSE calculates the RMSE of x and xa from dam_D objects.
        %
        %  - Input variable(s) -
        %  OBJ: A dam_D object        
        %
        %  NRS: The variable numbers that need to be taken into account. 
        %  For example:
        %       NRS=[1:3,6] calculates RMSE of variable numbers 1 to 3 and 6.        
        %
        %  - Output variable(s) -
        %  RMSEVAL: column vector containing the RMSE values of x and xa 
        %  as specified in NRS.        
        %
        %  - Construction -     
        %  RMSEVAL=RMSE(OBJ,NRS) calculates the RMSE of x and xa as 
        %  specified in NRS from dam_D objects.
        
        narginchk(2,2);  
            
        %error checking
            if isempty(obj)            
                error('DA:DataAssimilationModels:dam_D:dam_D:objEmpty','Can not calculate RMSE of an empty simulator object.')
            end
 
            if isempty(obj.x)
                error('DA:DataAssimilationModels:dam_D:dam_D:xEmpty','Can not calculate RMSE since x is empty.')
            end                
            
            % make nrs a row vector
            if size(nrs,1)>1; nrs=nrs';end;
            
            if min(nrs)<1 || (max(nrs)>size(obj.xa,1)) || (max(nrs)>size(obj.x,1))
                error('DA:DataAssimilationModels:dam_D:dam_D:nrsErr','RMSE: Index nrs exceeds object dimensions')
            end                        
        
            rmseval = sqrt( sum( (obj.x(nrs,:)-obj.xa(nrs,:)).^2,2) / size(obj.x,2) );
            
        end
        
        function mseval=mse(obj,nrs)
        %RMSE calculates the MSE of x and xa from dam_D objects.
        %
        %  - Input variable(s) -
        %  OBJ: A dam_D object        
        %
        %  NRS: The variable numbers that need to be taken into account. 
        %  For example:
        %       NRS=[1:3,6] calculates MSE of variable numbers 1 to 3 and 6.        
        %
        %  - Output variable(s) -
        %  MSEVAL: column vector containing the MSE values of x and xa 
        %  as specified in NRS.        
        %
        %  - Construction -     
        %  MSEVAL=MSE(OBJ,NRS) calculates the RMSE of x and xa as 
        %  specified in NRS from dam_D objects.
        
        narginchk(2,2);  
            
        %error checking
            if isempty(obj)            
                error('DA:DataAssimilationModels:dam_D:dam_D:objEmpty','Can not calculate MSE of an empty simulator object.')
            end
 
            if isempty(obj.x)
                error('DA:DataAssimilationModels:dam_D:dam_D:xEmpty','Can not calculate MSE since x is empty.')
            end                
            
            % make nrs a row vector
            if size(nrs,1)>1; nrs=nrs';end;
            
            if min(nrs)<1 || (max(nrs)>size(obj.xa,1)) || (max(nrs)>size(obj.x,1))
                error('DA:DataAssimilationModels:dam_D:dam_D:nrsErr','MSE: Index nrs exceeds object dimensions')
            end                        
        
            mseval =  sum( (obj.x(nrs,:)-obj.xa(nrs,:)).^2,2) / size(obj.x,2);
            
        end        
        
        function plot(obj,conf,nrs,varargin)
        %PLOT Overloaded plot method to plot dam_D objects.
        %
        %  - Input variable(s) -
        %  OBJ: A dam_D object 
        %
        %  CONF: String that defines the variable(s) to plot, the desired 
        %  scale of the x axis and the type of plot. 
        %  Can be 'u * **', 'x * **', 'xa * **', 'xf * **', 'Pa * **', 
        %  'Pf * **', 'y * **', 'xa-x * **', 'xaPa *', 'xaxPa *', 'xfPf *', 
        %  'xy *', 'xax *' and 'xaxy *' where * is either
        %   - 'k' to plot in function of the step nr index (default) or 
        %   - 't' to plot in function of the time.
        %  Likewise, a number can be added to define the kind of desired plot
        %  by defining **, where ** is either
        %   - 1 to create one plot of all variables or 
        %   - 2 to create a subplot per variable (2=default).
        %  Note1: 'Pa' & 'Pf' only plots the variances and 'xaPa' and 
        %  'xfPf' plot xa resprectively xf with their 95% conf. interval.
        %  Note2: ** can not be chosen for 'xy', 'xaPa', 'xaxPa *', 'xfPf'
        %  , 'xaxy' and 'xax' since they are only available with subplots.
        %  Note3: 'xa-x' plots the difference between est. and true state.
        %
        %  NRS: The variable numbers that need to be plotted. Maximum
        %  15 variables can be plotted at the same time.
        %  For example:
        %       NRS=[1:3,6] plots variable numbers 1 to 3 and 6.
        %
        %  VARARGIN: additional arguments can be added as in regular plot
        %  function. For example: plot(obj,'y',[5:9],'linewidth',3) plots
        %  the measuremnts 5 to 9 in linewidth 3.
        %
        %  - Output variable(s) -
        %  No
        %  
        %  - Construction -         
        %  TO PLOT ONE VARIABLE OR VARS WITH IDENTICAL NRS:
        %  PLOT(OBJ,CONF,NRS,VARARGIN) plots the requested variable CONF and
        %  more specifically the variable numbers provided in NRS while
        %  using VARARGIN to set the plot style.
        %  
        %  PLOT(OBJ,CONF,NRS) plots the requested variable CONF and
        %  more specifically the variable numbers provided in NRS. The plot
        %  style, legend and title of the plot are predefined.
        %
        %  TO PLOT A COMPARISON OF VARIABLES WITH DIFFERENT NRS: 'xy', 'xaxy'
        %  PLOT(OBJ,CONF,NRS1,NRS2,VARARGIN) plots the requested variable CONF
        %  and more specifically the variable numbers provided in NRS1 for x
        %  and NRS2 for y while using VARARGIN to set the plot style. 
        %  Both NRS have to contain an equal amount of variable numbers. 
        %  Plot style 1 is not supported.
        %  
        %  PLOT(OBJ,CONF,NRS1,NRS2) plots the requested variable CONF and
        %  more specifically the variable numbers provided in NRS1 for x
        %  and NRS2 for y. The plot style, legend and title of the plot  
        %  are predefined. Both NRS have to contain an equal amount of 
        %  variable numbers. Plot style 1 is not supported.            

            narginchk(3,inf);            
            ni=nargin;
            
                %Resolve outConf in bits
                remain = conf;
                plotxa=0;plotxf=0;plotx=0;plotu=0;ploty=0;plott=0;plotstyle=0;
                plotPa=0;plotPf=0;plotxaPa=0;plotxfPf=0;plotxaxPa=0;plotxaxy=0;
                plotxax=0;plotxy=0;plotxa_x=0;
                while true
                   [str, remain] = strtok(remain, ' '); %#ok<STTOK>
                   if isempty(str),  break;  end

                   if ~plotxa;plotxa=strcmp(str,'xa');end
                   if ~plotxf;plotxf=strcmp(str,'xf');end
                   if ~plotx;plotx=strcmp(str,'x');end
                   if ~plotu;plotu=strcmp(str,'u');end
                   if ~ploty;ploty=strcmp(str,'y');end
                   if ~plott;plott=strcmp(str,'t');end
                   if ~plotPa;plotPa=strcmp(str,'Pa');end                   
                   if ~plotPf;plotPf=strcmp(str,'Pf');end                   
                   if ~plotxaPa;plotxaPa=strcmp(str,'xaPa');end     
                   if ~plotxfPf;plotxfPf=strcmp(str,'xfPf');end                        
                   if ~plotxax;plotxax=strcmp(str,'xax');end   
                   if ~plotxy;plotxy=strcmp(str,'xy');end      
                   if ~plotxa_x;plotxa_x=strcmp(str,'xa-x');end  
                   if ~plotxaxPa;plotxaxPa=strcmp(str,'xaxPa');end  
                   if ~plotxaxy;plotxaxy=strcmp(str,'xaxy');end  
                   if plotstyle==0;plotstyle=strcmp(str,'1')+2*strcmp(str,'2');end
                end
            
            if ~plotxa && ~plotxf && ~plotx && ~plotu && ~ploty && ~plott && ~plotPa && ~plotPf && ~plotxaPa ...
               && ~plotxfPf  && ~plotxax  && ~plotxy  && ~plotxa_x && ~plotxaxPa && ~plotxaxy     
                error('DA:DataAssimilationModels:dam_D:dam_D:confErr','Please specify what you want to plot in ''conf''.')
            end
            
            % enforce plotstyle 2 when default, xy, xaPa, xfPf, xax
            if plotstyle==0 || plotxy || plotxaPa || plotxaxPa || plotxfPf || plotxax || plotxaxy
                plotstyle=2;
            end

            %handle two diff. vars
            if plotxy || plotxaxy
                nrs2 = varargin{1};
                varargin = varargin(2:end);
                ni=ni-1; %adjust nr of input arguments for transparancy
            else
                nrs2=nrs;
                if ni>3 
                    if isa(varargin{1},'double')           
                    error('DA:DataAssimilationModels:dam_D:dam_D:inputArgs','Two NRS is only allowed when plotting x and y.')
                    end
                end
            end
                        
            %error checking
            if isempty(obj)            
                error('DA:DataAssimilationModels:dam_D:dam_D:objEmpty','Can not plot empty simulator object.')
            end
            
            if isempty(obj.u) && plotu
                error('DA:DataAssimilationModels:dam_D:dam_D:uEmpty','Can not plot empty u.')
            end
         
            if isempty(obj.x) && (plotx || plotxy || plotxax || plotxaxy || plotxa_x || plotxaxPa)
                error('DA:DataAssimilationModels:dam_D:dam_D:xEmpty','Can not plot empty x.')
            end            

            if isempty(obj.xf) && (plotxf || plotxfPf)
                error('DA:DataAssimilationModels:dam_D:dam_D:xEmpty','Can not plot empty xf.')
            end         
            
            if isempty(obj.Pf) && (plotPf || plotxfPf)
                error('DA:DataAssimilationModels:dam_D:dam_D:xEmpty','Can not plot empty Pf.')
            end          
            
            if isempty(obj.Pa) && (plotPa || plotxaPa || plotxaxPa)
                error('DA:DataAssimilationModels:dam_D:dam_D:xEmpty','Can not plot empty Pf.')
            end                
            
            if isempty(obj.y) && (ploty || plotxy || plotxaxy)
                error('DA:DataAssimilationModels:dam_D:dam_D:xEmpty','Can not plot empty y.')
            end         
            
            if length(nrs)>15 || length(nrs2)>15
                error('DA:DataAssimilationModels:dam_D:dam_D:maxNrs','Can not plot more than 15 variables.')
            end                  
            
            if length(nrs) ~= length(nrs2)
                error('DA:DataAssimilationModels:dam_D:dam_D:Nrs1Nrs2','NRS1 and NRS2 need to define an equal amount of variables')
            end         
            
            % make nrs a row vector
            if size(nrs,1)>1; nrs=nrs';end;
            if size(nrs2,1)>1; nrs2=nrs2';end;
            
            if min(nrs)<1 || min(nrs2)<1 ...
                          || (max(nrs)>size(obj.xf,1) && (plotxf || plotxfPf) ) ...
                          || (max(nrs)>size(obj.xa,1) && (plotxa || plotxaPa || plotxaxPa || plotxax || plotxa_x || plotxaxy) ) ...
                          || (max(nrs)>size(obj.Pa,1) && (plotPa || plotxaPa || plotxaxPa) ) ...
                          || (max(nrs)>size(obj.Pf,1) && (plotPf || plotxfPf) ) ...                          
                          || (max(nrs)>size(obj.x,1) && (plotx || plotxy || plotxax || plotxa_x || plotxaxPa || plotxaxy) ) ...
                          || (max(nrs)>size(obj.y,1) && (ploty ) )...
                          || (max(nrs2)>size(obj.y,1) && (plotxy  || plotxaxy ) )...                          
                          || (max(nrs)>size(obj.u,1) && plotu) 
                error('DA:DataAssimilationModels:dam_D:dam_D:nrsErr','Index nrs exceeds object dimensions')
            end                
            
            maxPlotPnts=1000;
            if length(obj.k)>maxPlotPnts
                plotskipfact = floor(length(obj.k)/maxPlotPnts);
                dataind = 1:plotskipfact:length(obj.k);
                kind=obj.k(1):plotskipfact:obj.k(end);
            else
                kind = obj.k;
                dataind = 1:1:length(obj.k);
            end
            
            % define x axis scale and label
            if plott
                index = kind*obj.Ts;
                xLblStr = sprintf('%s',obj.TimeUnit);
            else
                index = kind;
                xLblStr = sprintf('k');
            end
            
            % define y axis label, title and legend text
            if plotxy;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - true states (x) & measurements (y)');end;
            if plotx;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - true states (x)');end;
            if plotxa || plotxaPa;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - analysis states (xa)');end;   
            if plotxaxPa;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - analysis states (xa) & true states (x)');end;                
            if plotxf || plotxfPf;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - forecast states (xf)');end;               
            if ploty;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - measurements (y)');end;
            if plotu;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - inputs (u)');end;
            if plotxax;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - analysis states (xa) & true states (x)');end;    
            if plotxa_x;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - state error (xa-x)');end;    
            if plotPf;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - forecast variances (Pf)');end;    
            if plotPa;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - analysis variances (Pa)');end;        
            if plotxaxy;yLblStr=' ';titleStr=strcat(obj.alg,'- Assimilation - analysis (xa), true states (x) & measurements (y)');end;            
            
            legendTxt=strtrim(cellstr(num2str((nrs)')));
            for i=1:length(legendTxt)
            	legendTxt{i}=char(yLblStr',legendTxt{i})';
            end    
            
            switch plotstyle
                case 1  %all in one plot
                    if ni==3
                        figure('Position', [300, 200, 800, 400]);hold on;  

                        if plotx
                            plot(index,obj.x(nrs,dataind),'linewidth',1);
                        elseif plotxf
                            plot(index,obj.xf(nrs,dataind),'linewidth',1);   
                        elseif plotxa
                            plot(index,obj.xa(nrs,dataind),'linewidth',1);                              
                        elseif ploty
                            plot(index,obj.y(nrs,dataind),'linewidth',1);               
                        elseif plotu
                            plot(index,obj.u(nrs,dataind),'linewidth',1); 
                        elseif plotPa
                            plot(index,getVars(obj.Pa,obj.covFull,nrs,dataind),'linewidth',1);     
                        elseif plotPf
                            plot(index,getVars(obj.Pf,obj.covFull,nrs,dataind),'linewidth',1);           
                        elseif plotxa_x
                            plot(index,obj.xa(nrs,dataind)-obj.x(nrs,dataind),'linewidth',1); 
                        end
                        axis tight;
                        axis 'auto y';
                        ylabel(yLblStr,'FontSize', 12);
                        title(titleStr,'FontSize', 13);                           
                        xlabel(xLblStr,'FontSize', 12); hold off; 
                        legend(legendTxt);  
                    else
                        if plotx   
                            plot(index,obj.x(nrs,dataind),varargin{:});
                        elseif plotxf   
                            plot(index,obj.xf(nrs,dataind),varargin{:});
                        elseif plotxa   
                            plot(index,obj.xa(nrs,dataind),varargin{:});                            
                        elseif ploty
                            plot(index,obj.y(nrs,dataind),varargin{:});                            
                        elseif plotu
                            plot(index,obj.u(nrs,dataind),varargin{:});   
                        elseif plotPa
                            plot(index,getVars(obj.Pa,obj.covFull,nrs,dataind),varargin{:});     
                        elseif plotPf
                            plot(index,getVars(obj.Pf,obj.covFull,nrs,dataind),varargin{:});           
                        elseif plotxa_x
                            plot(index,obj.xa(nrs,dataind)-obj.x(nrs,dataind),varargin{:});                             
                        end
                    end
                case 2  %subplots
                   
                    %dimensions of subplots made easy
                	[rows,cols]=numSubplots(length(nrs));                    
                    
                    if ni==3
                        figure('Position', [300, 200, 800, 400]);hold on; 

                        for i=1:length(nrs)
                            subplot(rows,cols,i);
                            if plotx
                                plot(index,obj.x(nrs(i),dataind),'b','linewidth',1);
                                titleStr=sprintf('x%s',num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);
                            elseif plotxa
                                plot(index,obj.xa(nrs(i),dataind),'b','linewidth',1);
                                titleStr=sprintf('xa%s',num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);                                   
                            elseif plotxa_x
                                plot(index,obj.xa(nrs(i),dataind)-obj.x(nrs(i),dataind),'b','linewidth',1);
                                titleStr=sprintf('xa%s -x%s',num2str(nrs(i)),num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);                                   
                            elseif plotxf
                                plot(index,obj.xf(nrs(i),dataind),'b','linewidth',1);
                                titleStr=sprintf('xf%s',num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);    
                            elseif plotPa
                                plot(index,getVars(obj.Pa,obj.covFull,nrs(i),dataind),'k','linewidth',1);
                                titleStr=sprintf('Pa%s',num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);       
                            elseif plotPf
                                plot(index,getVars(obj.Pf,obj.covFull,nrs(i),dataind),'k','linewidth',1);
                                titleStr=sprintf('Pf%s',num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);                                 
                            elseif ploty
                                plot(index,obj.y(nrs(i),dataind),'r','linewidth',1);    
                                titleStr=sprintf('y%s',num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);
                            elseif plotu
                                plot(index,obj.u(nrs(i),dataind),'g','linewidth',1); 
                                titleStr=sprintf('u%s',num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);
                            elseif plotxy
                                plot(index,obj.x(nrs(i),dataind),'b','linewidth',1);hold on;
                                plot(index,obj.y(nrs2(i),dataind),'r--','linewidth',1);    
                                titleStr=sprintf('x%s (blue) vs y%s (red)',num2str(nrs(i)),num2str(nrs2(i)));
                                title(titleStr,'FontSize', 10);
                            elseif plotxax
                                plot(index,obj.xa(nrs(i),dataind),'b','linewidth',1);hold on;
                                plot(index,obj.x(nrs(i),dataind),'g--','linewidth',1);    
                                titleStr=sprintf('xa%s (blue) vs x%s (green)',num2str(nrs(i)),num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);     
                            elseif plotxaxy
                                plot(index,obj.xa(nrs(i),dataind),'b','linewidth',1);hold on;
                                plot(index,obj.x(nrs(i),dataind),'g--','linewidth',1);    
                                plot(index,obj.y(nrs2(i),dataind),'r--','linewidth',1);                                    
                                titleStr=sprintf('xa%s (blue) - x%s (green) - y%s (red)',num2str(nrs(i)),num2str(nrs(i)),num2str(nrs2(i)));
                                title(titleStr,'FontSize', 10);                                  
                            elseif plotxaPa
                                PaCI=1.96*(sqrt(getVars(obj.Pa,obj.covFull,nrs(i),dataind)));
                                plot(index,obj.xa(nrs(i),dataind),'b','linewidth',1);hold on;
                                plot(index,obj.xa(nrs(i),dataind)-PaCI,'m--','linewidth',1);    
                                plot(index,obj.xa(nrs(i),dataind)+PaCI,'m--','linewidth',1);    
                                titleStr=sprintf('xa%s (blue)',num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);    
                            elseif plotxfPf
                                PfCI=1.96*(sqrt(getVars(obj.Pf,obj.covFull,nrs(i),dataind)));
                                plot(index,obj.xf(nrs(i),dataind),'b','linewidth',1);hold on;
                                plot(index,obj.xf(nrs(i),dataind)-PfCI,'m--','linewidth',1);    
                                plot(index,obj.xf(nrs(i),dataind)+PfCI,'m--','linewidth',1);    
                                titleStr=sprintf('xf%s (blue)',num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);         
                            elseif plotxaxPa
                                PaCI=1.96*(sqrt(getVars(obj.Pa,obj.covFull,nrs(i),dataind)));
                                plot(index,obj.xa(nrs(i),dataind),'b','linewidth',1);hold on;
                                plot(index,obj.x(nrs(i),dataind),'g--','linewidth',1);                        
                                plot(index,obj.xa(nrs(i),dataind)-PaCI,'m--','linewidth',1);    
                                plot(index,obj.xa(nrs(i),dataind)+PaCI,'m--','linewidth',1);    
                                titleStr=sprintf('xa%s (blue) - x%s (green)',num2str(nrs(i)),num2str(nrs(i)));
                                title(titleStr,'FontSize', 10);                                   
                            end
                            
                            axis tight;
                            axis 'auto y';                            
                            ylabel(yLblStr,'FontSize', 10);  
                            xlabel(xLblStr,'FontSize', 10);                          
                        end        
                    else
                        for i=1:length(nrs)
                            subplot(rows,cols,i);
                            if plotx
                                plot(index,obj.x(nrs(i),dataind),varargin{:});
                            elseif plotxa
                                plot(index,obj.xa(nrs(i),dataind),varargin{:});                      
                            elseif plotxa_x
                                plot(index,obj.xa(nrs(i),dataind)-obj.x(nrs(i),:),varargin{:});                           
                            elseif plotxf
                                plot(index,obj.xf(nrs(i),dataind),varargin{:}); 
                            elseif plotPa
                                plot(index,getVars(obj.Pa,obj.covFull,nrs(i),dataind),varargin{:});  
                            elseif plotPf
                                plot(index,getVars(obj.Pf,obj.covFull,nrs(i),dataind),varargin{:});                            
                            elseif ploty
                                plot(index,obj.y(nrs(i),dataind),varargin{:});
                            elseif plotu
                                plot(index,obj.u(nrs(i),dataind),varargin{:});
                            elseif plotxy
                                plot(index,obj.x(nrs(i),dataind),'b',varargin{:});hold on;
                                plot(index,obj.y(nrs2(i),dataind),'r--',varargin{:});   
                            elseif plotxax
                                plot(index,obj.xa(nrs(i),dataind),'b',varargin{:});hold on;
                                plot(index,obj.x(nrs(i),dataind),'g--',varargin{:});     
                            elseif plotxaxy
                                plot(index,obj.xa(nrs(i),dataind),'b',varargin{:});hold on;
                                plot(index,obj.x(nrs(i),dataind),'g--',varargin{:});    
                                plot(index,obj.y(nrs2(i),dataind),'r--',varargin{:});                                                              
                            elseif plotxaPa
                                PaCI=1.96*(sqrt(getVars(obj.Pa,obj.covFull,nrs(i),dataind)));
                                plot(index,obj.xa(nrs(i),dataind),'b',varargin{:});hold on;
                                plot(index,obj.xa(nrs(i),dataind)-PaCI,'m--',varargin{:});  
                                plot(index,obj.xa(nrs(i),dataind)+PaCI,'m--',varargin{:});     
                            elseif plotxfPf
                                PfCI=1.96*(sqrt(getVars(obj.Pf,obj.covFull,nrs(i),dataind)));
                                plot(index,obj.xf(nrs(i),dataind),'b',varargin{:});hold on;
                                plot(index,obj.xf(nrs(i),dataind)-PfCI,'m--',varargin{:});  
                                plot(index,obj.xf(nrs(i),dataind)+PfCI,'m--',varargin{:});          
                            elseif plotxaxPa
                                PaCI=1.96*(sqrt(getVars(obj.Pa,obj.covFull,nrs(i),dataind)));
                                plot(index,obj.xa(nrs(i),dataind),'b',varargin{:});  hold on;
                                plot(index,obj.x(nrs(i),dataind),'g--',varargin{:});                        
                                plot(index,obj.xa(nrs(i),dataind)-PaCI,'m--',varargin{:});    
                                plot(index,obj.xa(nrs(i),dataind)+PaCI,'m--',varargin{:});                                   
                            end                
                        end    
                    end
            end        
        
        end
        
        
    end %methods
    
    methods (Access = private)    
        
        function value=checkConsist(obj)
        %CHECKCONSIST Checks if object is consistent.
            value=0;
            if ~isempty(obj)

                if (~isempty(obj.Pf) || ~isempty(obj.Pa)) && ~isequal(obj.covFull,0) && ~isequal(obj.covFull,1)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','covFull must be set to one or zero.')                     
                end
                
                if ~isempty(obj.xf) && size(obj.xf,1)~=size(obj.xa,1)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','xf and xa must have same number of rows.') 
                end
                
                if ~isempty(obj.Pf) && size(obj.Pf,1)~=size(obj.xa,1)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','Pf and xa must have same number of rows.') 
                end                
                
                if ~isempty(obj.Pa) && size(obj.Pa,1)~=size(obj.xa,1)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','Pa and xa must have same number of rows.') 
                end                           
                
                if ~isempty(obj.xf) && size(obj.xf,2)~=size(obj.xa,2)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','xf and xa must have same number of columns.') 
                end

                if ~isempty(obj.u) && size(obj.u,2)~=size(obj.xa,2)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','u and xa must have same number of columns.') 
                end                
                
                if ~isempty(obj.y) && size(obj.y,2)~=size(obj.xa,2)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','y and xa must have same number of columns.') 
                end                       
                
                if ~isempty(obj.Pf) && ~obj.covFull && size(obj.Pf,2)~=size(obj.xa,2)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','Pf and xa must have same number of columns.') 
                end       
                
                if ~isempty(obj.Pf) && obj.covFull && size(obj.Pf,3)~=size(obj.xa,2)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','Pf and xa dimensions do not match.') 
                end                  
          
                if ~isempty(obj.Pa) && ~obj.covFull && size(obj.Pa,2)~=size(obj.xa,2)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','Pa and xa must have same number of columns.') 
                end       
                
                if ~isempty(obj.Pa) && obj.covFull && size(obj.Pa,3)~=size(obj.xa,2)
                   error('DA:DataAssimilationModels:dam_D:dam_D:dimMismatch','Pa and xa dimensions do not match.') 
                end                                      
                
                value=1;
            end
        
        end
    
    end %methods    
    
end