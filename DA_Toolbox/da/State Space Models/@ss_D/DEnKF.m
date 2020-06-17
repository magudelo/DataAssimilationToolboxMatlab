function daModel=DEnKF(ssObj,varargin)
%DEnKF performs Deterministic Ensemble Kalman Filter algorithm on a 
%    Discrete Linear State-Space object
%
%  - Input variable(s) -
%  SSOBJ: discrete State-Space object. (type 'help ss_D')
%
%  Y: matrix of measurements where every column is a measurement vector.
%  The first measurement has to correspond to x0, the initial state at step
%  k0 of ssObj.
%
%  U: matrix of inputs where every column is a input vector. Can be left 
%  empty ([ ]) if desired. 
%
%  SIMMODEL: a discrete time State-Space simulation model of type sim_D.
%  The matrices Y and U are obtained from this model.
%
%  OUTCONF: string that indicates which variables must reside in the
%  data assimilation model. Possible string values contain 'xf', 'Pa', 
%  'Pf', 'y', 'u', 'cov', 'var'. 'xa' is always returned.
%  For example: 'xf y' retains xa, xf and y matrices.
%  Default: 'xa'
%
%  OPTIONSCELL: cell that contains additional options. The Deterministic 
%  Ensemble Kalman  Filter has the option NENSEMBLES or INITENSEMBLES. 
%  NENSEMBLES indicates the number of ensembles to use. When NENSEMBLES is 
%  used, the ensembles  are selected randomly from the X0 noise model of  
%  the SSOBJ. INITENSEMBLES contains user defined initial ensembles.
%  When neither NENSEMBLES or INITENSEMBLES are configured NENSEMBLES is
%  set equal to the amount of states of SSOBJ with a minimum of 5 and a 
%  maximum of 50.
%
%  - Output variable(s) -
%  DAMODEL: a discrete time State-Space data assimilation model of type dam_D.
%  
%  - Construction -          
%  DAMODEL=DEnKF(SSOBJ,Y,U,OUTCONF,OPTIONSCELL) performs the  
%  Ensemble Kalman   Filter algorithm on a Discrete State-Space object.
%  The parameters U, OUTCONF and OPTIONSCELL can be omitted if desired but
%  the order of the parameters must remain.
%  For example: DEnKF(SSOBJ,Y,U,OPTIONSCELL) and DEnKF(SSOBJ,Y,OUTCONF) are
%  allowed, but DEnKF(SSOBJ,OUTCONF,Y) is not allowed
%
%  DAMODEL=DEnKF(SSOBJ,SIMMODEL,OUTCONF,OPTIONSCELL) performs the   
%  Deterministic Ensemble Kalman Filter algorithm on a Discrete State-Space   
%  object using data from the simulation model SIMMODEL.
%  The parameters OUTCONF and OPTIONSCELL can be omitted if desired but
%  the order of the parameters must remain.
%  For example: DEnKF(SSOBJ,SIMMODEL,OPTIONSCELL) and DEnKF(SSOBJ,Y,OUTCONF) 
%  are allowed, but DEnKF(SSOBJ,SIMMODEL,OPTIONSCELL,OUTCONF) is not allowed
%
% References: 
% A deterministic formulation of the ensemble Kalman filter: 
% an alternative to ensemble square root filters. PAVEL SAKOV,PETER R. OKE
% Tellus A, Volume 60, Issue 2, pages 361–371, March 2008

    %======================================%
    %   RESOLVE ARGS & INITIAL CHECKINGS   %
    %======================================%

    % check amount of input arguments
	narginchk(2,6);
	ni = nargin; 
    
    %%GENERAL call 'da(ssObj,alg,varargin)': second argument is a cell       
	if ni==2 && isa(varargin{1},'cell')
        varargin=varargin{1};
        if length(size(varargin))> 2 || (size(varargin,1)>1 && size(varargin,2)>1) 
            error('DA:StateSpaceModels:ss_D:DEnKF:cellDim','Cell must be vector sized.')
        else
            ni=1+max(size(varargin,1),size(varargin,2));
        end
	end    
    
    % Retrieve k0 from ssObj       
    k0=ssObj.k0;
    
    % Deterministic Ensemble Kalman Filter should be applied on additive noise models
    if isa(ssObj,'ss_DNL') || isa(ssObj,'ss_DNL_AMN')
        warning('DA:StateSpaceModels:ss_D:DEnKF:ssWarning','Deterministic Ensemble Kalman Filter should only be applied on additive noise models.')
    end    
    
    % Deterministic Ensemble Kalman Filter should be applied on zero mean gaussian noise       
    if ~iszeromean(mean(ssObj.w,k0)) || ~iszeromean(mean(ssObj.v,k0)) || ~isa(ssObj.w,'nm_gauss') || ~isa(ssObj.v,'nm_gauss')
        warning('DA:StateSpaceModels:ss_D:DEnKF:nonZeroMean','Noise models should be zero mean and Gaussian for DEnKF.')
    end
    
    % Simulation object is provided: save simObj, shift cell content and 
    % place simObj.y and simObj.u in first and second place of varargin
    xtrue=[];simTrue=0;
    if isa(varargin{1},'sim_D')
        simObj=varargin{1};                 % first object is simObj
        varargin(end+1)={0};                % make last el. zero
        varargin=circshift(varargin,[0 1]); % shift elements to the right        
        varargin{1}=simObj.y;               
        varargin{2}=simObj.u;
        simTrue=1;
        ni=ni+1;                            % ni is increased
    end

    %Initial settings     
    optionsCell=[];
    outConf='xa';
    u=[];
    
    %Resolve arguments    
    if isa(varargin{1},'double')
        y=varargin{1};                                  %First argument must be y  
        if ni>2
            if isa(varargin{2},'double')
                u=varargin{2};                          %Second argument if double is u 
                if ni > 3 
                    if isa(varargin{3},'char')          %Third argument if char is outConf 
                        outConf=varargin{3};
                        if ni>4
                            if isa(varargin{4},'cell')  %Fourth argument if cell is optionsCell
                                optionsCell=varargin{4};
                            end
                        end
                    elseif isa(varargin{3},'cell')      %Third argument if not outConf is optionsCell  
                        optionsCell=varargin{4};
                    end
                end 
            elseif isa(varargin{2},'char')              %Second argument if not u is outConf  
                outConf=varargin{2};    
                if ni>3
                    if isa(varargin{3},'cell')          %Third argument if not u is optionsCell
                        optionsCell=varargin{4};
                    end
                end                
            elseif isa(varargin{2},'cell')              %Second argument if not u and not outConf is optionsCell  
                optionsCell=varargin{4};
            end
        end
    else
        error('DA:StateSpaceModels:ss_D:DEnKF:classMismatch','Class mismatch for y: expected sim_D object or double.')
    end

    %Resolve optionsCell
    if ~isempty(optionsCell)
        ni=max([size(optionsCell,1) size(optionsCell,2)]);
        if ni>0
            n=size(optionsCell{1});
            if isa(optionsCell{1},'double')&&length(n)==2&&n(1)==1&&n(2)==1         %scalar
                nEnsembles=optionsCell{1};
                initEnsembles=0;
            elseif isa(optionsCell{1},'double')&&length(n)==2&&n(1)>1&&n(2)>1       %matrix 
                initEnsembles=optionsCell{1};
                nEnsembles = size(initEnsembles,2);
                if size(initEnsembles,1)~=ssObj.x0
                    error('DA:StateSpaceModels:ss_D:DEnKF:initEnsembles','Number of states of initEnsembles is incorrect.')
                end  
            elseif isempty(optionsCell{1})
                nEnsembles = size(ssObj.x0,1);
                if nEnsembles>50;nEnsembles=50;end
                if nEnsembles<5;nEnsembles=5;end      
                initEnsembles=0;
            else
                error('DA:StateSpaceModels:ss_D:DEnKF:argMismatch','Optionscell must contain nEnsembles or initEnsembles as first argument.')
            end 
        else
            %Nr of ensembles is taken equal to amount of states (min5,max50)
            nEnsembles = size(ssObj.x0,1);
            if nEnsembles>50;nEnsembles=50;end
            if nEnsembles<5;nEnsembles=5;end   
            initEnsembles=0;
        end
    else
        %Nr of ensembles is taken equal to amount of states (min5,max50)
        nEnsembles = size(ssObj.x0,1);
        if nEnsembles>50;nEnsembles=50;end
        if nEnsembles<5;nEnsembles=5;end   
        initEnsembles=0;
    end
    if nEnsembles<2
        error('DA:StateSpaceModels:ss_D:DEnKF:nEnsembles','Number of ensembles must be minimal equal to two.')
    end    
    
    %Resolve outConf in bits
    remain = outConf;
    xtrueRet=0;xfRet=0;PfRet=0;PaRet=0;yRet=0;uRet=0;covFull=0;
    while true
       [str, remain] = strtok(remain, ' '); %#ok<STTOK>
       if isempty(str),  break;  end
       
       if ~xtrueRet;xtrueRet=strcmp(str,'x');end
       if ~xfRet;xfRet=strcmp(str,'xf');end
       if ~PfRet;PfRet=strcmp(str,'Pf');end
       if ~PaRet;PaRet=strcmp(str,'Pa');end
       if ~yRet;yRet=strcmp(str,'y');end
       if ~uRet;uRet=strcmp(str,'u');end      
       if ~covFull;covFull=strcmp(str,'cov') && ~strcmp(str,'var');end      
    end
    
    %Determine amount of samples based on provided measurements and inputs
    %(smallest amount counts)    
    samples=size(y,2);
    
    if ~isequal(u,0) && ~isempty(u) && size(u,2)~=size(y,2)
        warning('DA:StateSpaceModels:ss_D:DEnKF:uSize','Amount of inputs u is not consistent with amount of measurements y. Smallest amount is used.')
        if size(u,2)<size(y,2);samples=size(u,2);end
        if size(u,2)>size(y,2);samples=size(y,2);end
    end        
    
    %Timings        
	Ts=ssObj.Ts;
	kIndex=k0:1:(k0+samples-1);	

    %======================================================%
    %    DETERMINISTIC ENSEMBLE KALMAN FILTER ALGORITHM    %
    %======================================================%
    
    %Initial ensemble: if empty create random ensemble    
    if initEnsembles==0 || isempty(initEnsembles)
        initEnsembles=sample(ssObj.x0,nEnsembles);
    end     
    
    %Initial estimates & pre-allocate memory
    xhat_ens=initEnsembles;
    xhat = mean(xhat_ens,2);
    nStates=size(xhat,1);
	nOutputs=size(y,1);
    yf_ens=zeros(nOutputs,nEnsembles);
    if PfRet || PaRet
        Ex=xhat_ens-repmat(xhat,1,nEnsembles);      %ensemble error matrix
        [Phat] = P_Save_Stats(Ex,covFull);
    end
    
    %pre-allocate memory of data arrays & fill in first sample
    if xfRet; xfDa = zeros(nStates,samples);xfDa(:,1)=xhat;end;
    if ~xfRet;xfDa=[];end;
    xaDa = zeros(nStates,samples);xaDa(:,1)=xhat;
    if PfRet && covFull; PfDa = zeros(nStates,nStates,samples);PfDa(:,:,1)=Phat; end;
    if PfRet && ~covFull; PfDa = zeros(nStates,samples);PfDa(:,1)=Phat; end;
    if ~PfRet;PfDa=[];end;
    if PaRet && covFull; PaDa = zeros(nStates,nStates,samples);PaDa(:,:,1)=Phat;end;
    if PaRet && ~covFull; PaDa = zeros(nStates,samples);PaDa(:,1)=Phat;end;
    if ~PaRet;PaDa=[];end;
    
    for i=2:samples
      
        %OBTAIN VALUES 
        muV=mean(ssObj.v,kIndex(i));          
        
        %FORECAST STEP
        for j=1:nEnsembles  %forecast every ensemble
            xhat_ens(:,j) = eval_ftot(ssObj,xhat_ens(:,j),kIndex(i-1),getU(u,i-1),[]);  % xf_ens (nxN)
            yf_ens(:,j) = eval_htot(ssObj,xhat_ens(:,j),kIndex(i),getU(u,i),muV);     % (mxN)
        end
        
        xhat=mean(xhat_ens,2);                        	%estimate xf: mean of ensembles
        yf=mean(yf_ens,2);                              %estimate yf: mean of meas. ensembles
                    
        Ex=xhat_ens-repmat(xhat,1,nEnsembles);          %ensemble perturbation matrix (nxN)
        Ey=yf_ens-repmat(yf,1,nEnsembles);              %ensemble error matrix of output error (mxN)
               
        Pxy=Ex*Ey' /(nEnsembles-1);                             %approximation of Pxy (n x m)
        Pyy=Ey*Ey' /(nEnsembles-1) + cov(ssObj.v,kIndex(i)); 	%approximation of Pyy=Ey-*Ey'/(N-1)+R (m x m) 

        if xfRet; xfDa(:,i)=xhat;end;                           %SAVE xf        
        if PfRet; [Phat] = P_Save_Stats(Ex,covFull);end;        %approximation of Pf (n x n)           
        if PfRet && covFull; PfDa(:,:,i)=Phat;end;              %SAVE Pf
        if PfRet && ~covFull; PfDa(:,i)=Phat;end;               %SAVE Pf           
        
        %KALMAN GAIN
        K = Pxy /Pyy;                                   % (n x m)   

        %ANALYSIS STEP
%         for j=1:nEnsembles              
%             Ex(:,j)=Ex(:,j)-0.5*K*eval_htot(ssObj,Ex(:,j),kIndex(i),getU(u,i),muV);                 %ensemble error matrix
%         end            
%         xhat = xhat + K*( y(:,i) - eval_htot(ssObj,xhat,kIndex(i),getU(u,i),muV));
%         xhat_ens = Ex + repmat(xhat,1,nEnsembles);
% 
%         if PaRet;  Phat = (Ex*Ex')/(nEnsembles-1);end;    %approximation of Pa (n x n)  
        
        % Faster:        
        xaTmp = 0.5*K*( y(:,i) - eval_htot(ssObj,xhat,kIndex(i),getU(u,i),muV));                                %0.5K(d-Hxf)
        xhat_ens = xhat_ens + 0.5*K*( repmat(y(:,i),1,nEnsembles) - yf_ens) + repmat(xaTmp,1,nEnsembles);       %xa_ens   
        xhat = xhat + 2*xaTmp;                                                                                  %xa

        if PaRet                                        
            Ex=xhat_ens-repmat(xhat,1,nEnsembles);          %ensemble error matrix               
            [Phat] = P_Save_Stats(Ex,covFull);       	%approximation of Pa (n x n) 
        end                            
        xaDa(:,i)=xhat;                                     %SAVE xa        
        if PaRet && covFull; PaDa(:,:,i)=Phat;end;          %SAVE Pa
        if PaRet && ~covFull; PaDa(:,i)=Phat;end;           %SAVE Pa    
        
    end
      
    %MAKE UNREQUIRED DATA EMPTY    
	if ~uRet;u=[];end;
    if ~yRet;y=[];end;
    if ~xtrueRet;xtrue=[];end;
        
    if ~isempty(u);u=u(:,1:samples);end;
    if ~isempty(y);y=y(:,1:samples);end;
    if xtrueRet && simTrue;xtrue=simObj.x(:,1:samples);end;    

    % + to change to double
    daModel = dam_D('DEnKF',xtrue,xfDa,xaDa,PfDa,PaDa,u ,y ,kIndex,+(covFull),Ts,ssObj.TimeUnit);  
                    
end

    
    
  