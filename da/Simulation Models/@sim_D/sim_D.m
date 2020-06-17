classdef sim_D
    %SIM_D Constructs a simulation model
    %
    %  SIM_D Properties:
    %
    %     x - Simulated states. Matrix of simulated states where each
    %     column is a state. Can be left empty ([ ]) if desired.
    %
    %     u - Inputs. Matrix of inputs corresponding to the simulated 
    %     states where each column is an input. Can be left empty ([ ]) if 
    %     desired.        
    %
    %     y - Simulated measurements. Matrix of simulated measurements 
    %     where each column is a measurement. (must be non-empty)   
    %
    %     k - step number array. Each element contains the step number
    %     corresponding to each column of x,u and y. For example: the
    %     second element of k corresponds to the step number of column two.
    %     (must be non-empty) 
    %
    %     Ts - sample time. (default is 1)   
    %
    %     TimeUnit - String representing the unit of the time variable. 
    %       TimeUnit can take the following values: 'nanoseconds',
    %       'microseconds', 'milliseconds', 'seconds', 'minutes', 'hours',
    %       'days', 'weeks', 'months', 'years'.
    %       Default: 'seconds'     
    %
    %     l - Total amount of iterations.  (read-only)
    %
    %  SIM_D Construction:
    %     OBJ = SIM_D(X,U,Y,K,TS,TIMEUNIT) creates an object OBJ
    %     representing a simulation object.
    %
    %     OBJ = SIM_D(X,U,Y,K,TS) creates an object OBJ representing a 
    %     simulation object with TIMEUNIT='seconds'.
    %
    %     OBJ = SIM_D(X,U,Y,K) creates an object OBJ representing a 
    %     simulation object with TIMEUNIT='seconds' and TS=1.
    %    
    %  SIM_D Methods:
    %     plot - Overloaded plot method to plot sim_D objects. 
    % 
    
	properties (SetAccess = private)

        x;              % Simulated states matrix
        u;              % Inputs matrix
        y;              % Simulated measurements meatrix(must be non-empty)  
        k;              % Step number array (must be non-empty)     
        Ts;             % Sample time
        TimeUnit;       % String representing the unit of the time variable.      
        
    end

	properties (Dependent = true, SetAccess = private)
        l;              % Total amount of iterations.
    end       
    
    methods
        
        function obj = sim_D(x,u,y,k,Ts,TimeUnit)
        % SIM_D Constructor 
            narginchk(4,6);
            
            ni=nargin;
            
            obj.x=x;            
            obj.u=u;
            obj.y=y;
            obj.k=k;           

            if ni<5;Ts=1;end;            
            if ni<6;TimeUnit='seconds';end;                
            if isempty(Ts);Ts=1;end; 
            if isempty(TimeUnit);TimeUnit='seconds';end;     
            
            obj.Ts=Ts;
            obj.TimeUnit=TimeUnit;
            
        end            
        
        function obj = set.x(obj,value)
            
            value=checkArgs(value,'x');     	%check x                    
            obj.x=value;
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
            
            value=checkArgs(value,'Ts');      	%check Ts             
            obj.Ts=value;
            checkConsist(obj);
            
        end    
        
        function obj = set.TimeUnit(obj,value)
            
            value=checkArgs(value,'TimeUnit');	%check TimeUnit                   
            obj.TimeUnit=value;
            checkConsist(obj);
            
        end                
        
        function value = get.l(obj)
            
            value = max([size(obj.x,2) size(obj.u,2) size(obj.y,2) length(obj.k)]);
            
        end                
        
        function display(obj)
        %DISPLAY Displays the object summary.
        %
        %  - Input variable(s) -
        %  OBJ: a sim_D object     
        %
        %  - Construction -
        %  DISPLAY(OBJ) displays the sim_D object summary.
        
        	disp(' ');
            disp([inputname(1),' = '])
            disp(' ');          
            if ~isempty(obj)              
                disp('Discrete Simulator Data Object.')
                disp('-------------------------------')
                disp('Properties:');
                if ~isempty(obj.x)
                    fprintf('States (.x) \t \t = \t[%i x %i %s] \n',size(obj.x,1),size(obj.x,2),class(obj.x));
                end
                if ~isempty(obj.u)                
                    fprintf('Inputs (.u) \t \t = \t[%i x %i %s] \n',size(obj.u,1),size(obj.u,2),class(obj.u));
                end
                fprintf('Measurements (.y) \t = \t[%i x %i %s] \n',size(obj.y,1),size(obj.y,2),class(obj.y));
                fprintf('Step nr index (.k) \t = \t[%i x %i %s] \n',size(obj.k,1),size(obj.k,2),class(obj.k));
                fprintf('Samples (.l) \t \t = \t%i \n',obj.l);                                
                fprintf('Sample Time(.Ts) \t = \t%i \n',obj.Ts);                
                fprintf('Time Unit(.TimeUnit) = \t%s \n',obj.TimeUnit);                
            else
                disp('Empty Discrete Simulator Data Object.')                
            end
        end   
        
        function value=isempty(obj)
        %ISEMPTY Check if object is empty.
        %
        %  - Input variable(s) -
        %  OBJ: a sim_D object   
        %
        %  - Output variable(s) -
        %  VALUE: indicates if object is empty (=1) or not empty (=0).
        %  
        %  - Construction -          
        %  VALUE = ISEMPTY(OBJ) checks if object is empty or not. The
        %  object is considered empty when y or k are empty.
        %  
        
            if (~isempty(obj.y) && ~isempty(obj.k))
                value=0;
            else
                value=1;
            end

        end          
        
        function plot(obj,conf,nrs,varargin)
        %PLOT Overloaded plot method to plot sim_D objects.
        %
        %  - Input variable(s) -
        %  OBJ: A sim_D object 
        %
        %  CONF: String that defines the variable(s) to plot, the desired 
        %  scale of the x axis and the type of plot. 
        %  Can be 'u * **', 'x * **', 'y * **', or 'xy *' where * is either
        %   - 'k' to plot in function of the step nr index (default) or 
        %   - 't' to plot in function of the time.
        %  Likewise, a number can be added to define the kind of desired plot
        %  by defining **, where ** is either
        %   - 1 to create one plot of all variables or 
        %   - 2 to create a subplot per variable (2=default).
        %  Note that ** can not be chosen for 'xy' since it only supports
        %  subplots.
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
        %  TO PLOT ONE VARIABLE:
        %  PLOT(OBJ,CONF,NRS,VARARGIN) plots the requested variable CONF and
        %  more specifically the variable numbers provided in NRS while
        %  using VARARGIN to set the plot style.
        %  
        %  PLOT(OBJ,CONF,NRS) plots the requested variable CONF and
        %  more specifically the variable numbers provided in NRS. The plot
        %  style, legend and title of the plot are predefined.
        %
        %  TO PLOT A COMPARISON OF TWO VARIABLES:
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
                plotxy=0;ploty=0;plotx=0;plotu=0;plott=0;plotstyle=0;
                while true
                   [str, remain] = strtok(remain, ' '); %#ok<STTOK>
                   if isempty(str),  break;  end

                   if ~plotxy;plotxy=strcmp(str,'xy');end
                   if ~ploty;ploty=strcmp(str,'y');end
                   if ~plotx;plotx=strcmp(str,'x');end
                   if ~plotu;plotu=strcmp(str,'u');end
                   if ~plott;plott=strcmp(str,'t');end
                   if plotstyle==0;plotstyle=strcmp(str,'1')+2*strcmp(str,'2');end
                end
            
            if ~plotxy && ~plotx && ~ploty && ~plotu && ~plott
                error('DA:SimulatorModels:sim_D:sim_D:confErr','Please specify what you want to plot in ''conf''.')
            end
            
            % enforce plotstyle 2 as default and xy
            if plotstyle==0;plotstyle=2;end;
            if plotxy;plotstyle=2;end;

            %handle two diff. vars
            if plotxy
                nrs2 = varargin{1};
                varargin = varargin(2:end);
                ni=ni-1; %adjust nr of input arguments for transparancy
            else
                nrs2=nrs;
                if ni>3 
                    if isa(varargin{1},'double')           
                    error('DA:SimulatorModels:sim_D:sim_D:inputArgs','Two NRS is only allowed when plotting x and y.')
                    end
                end
            end
                        
            %error checking
            if isempty(obj)            
                error('DA:SimulatorModels:sim_D:sim_D:objEmpty','Can not plot empty simulator object.')
            end
            
            if isempty(obj.u) && plotu
                error('DA:SimulatorModels:sim_D:sim_D:uEmpty','Can not plot empty u.')
            end
         
            if isempty(obj.x) && (plotx || plotxy)
                error('DA:SimulatorModels:sim_D:sim_D:xEmpty','Can not plot empty x.')
            end            
            
            if length(nrs)>15 || length(nrs2)>15
                error('DA:SimulatorModels:sim_D:sim_D:maxNrs','Can not plot more than 15 variables.')
            end                  
            
            if length(nrs) ~= length(nrs2)
                error('DA:SimulatorModels:sim_D:sim_D:Nrs1Nrs2','NRS1 and NRS2 need to define an equal amount of variables')
            end         
            
            % make nrs a row vector
            if size(nrs,1)>1; nrs=nrs';end;
            if size(nrs2,1)>1; nrs2=nrs2';end;
            
            if min(nrs)<1 || min(nrs2)<1 ...                     
                          || (max(nrs)>size(obj.x,1) && (plotx || plotxy) ) ...
                          || (max(nrs)>size(obj.y,1) && (ploty ) )...
                          || (max(nrs2)>size(obj.y,1) && (plotxy ) )...                          
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
            if plotxy;yLblStr=' ';titleStr='Simulator - true states x & measurements y';end;
            if plotx;yLblStr=' ';titleStr='Simulator - true states x';end;
            if ploty;yLblStr=' ';titleStr='Simulator - measurements y';end;
            if plotu;yLblStr=' ';titleStr='Simulator - inputs u';end;
            
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
                        elseif ploty
                            plot(index,obj.y(nrs,dataind),'linewidth',1);               
                        elseif plotu
                            plot(index,obj.u(nrs,dataind),'linewidth',1); 
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
                        elseif ploty
                            plot(index,obj.y(nrs,dataind),varargin{:});                            
                        elseif plotu
                            plot(index,obj.u(nrs,dataind),varargin{:});   
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
                            elseif ploty
                                plot(index,obj.y(nrs(i),dataind),varargin{:});  
                            elseif plotu
                                plot(index,obj.u(nrs(i),dataind),varargin{:});
                            elseif plotxy
                                plot(index,obj.x(nrs(i),dataind),'b',varargin{:});hold on;
                                plot(index,obj.y(nrs2(i),dataind),'r--',varargin{:});  
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
            
                if ~isempty(obj.x) && ( size(obj.x,2)~=size(obj.y,2) || size(obj.x,2) ~= size(obj.k,2) )
                   error('DA:SimulatorModels:sim_D:sim_D:dimMismatch','x,y and k must have same number of columns.') 
                end
                
                if ~isempty(obj.u) && (size(obj.u,2) ~= size(obj.k,2))
                   error('DA:SimulatorModels:sim_D:sim_D:dimMismatch','u and x,y,k must have same number of columns.') 
                end               
                
                value=1;
            end
        
        end
    
    end %methods    
    
end