function daModel=UKF(ssObj,varargin)
%UKF performs Unscented Kalman Filter algorithm on a Discrete Linear 
%    State-Space object
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
%  OPTIONSCELL: cell that contains additional options - {KAPPA,SPTWICE} or
%  {KAPPA} are possibilities where 
%  * KAPPA: used to reduce the higher-order errors of the mean and
%    covariance approximations. (Default=0) For example, if z is Gaussian  
%    then KAPPA = 3 - n will minimize some of the errors in the fourth-order 
%    terms in the mean and covariance approximation.
%  * SPTWICE: reuse Sigma Points of forecast step in analysis step (=0) 
%    or recalculate Sigma Points in analysis step (=1). Default=0. (saves
%    computional effort at the cost of performance)
%
%  - Output variable(s) -
%  DAMODEL: a discrete time State-Space data assimilation model of type dam_D.
%  
%  - Construction -          
%  DAMODEL=UKF(SSOBJ,Y,U,OUTCONF,OPTIONSCELL) performs the Unscented Kalman  
%  Filter algorithm on a Discrete State-Space object.
%  The parameters U, OUTCONF and OPTIONSCELL can be omitted if desired but
%  the order of the parameters must remain.
%  For example: UKF(SSOBJ,Y,U,OPTIONSCELL) and UKF(SSOBJ,Y,OUTCONF) are
%  allowed, but UKF(SSOBJ,OUTCONF,Y) is not allowed
%
%  DAMODEL=UKF(SSOBJ,SIMMODEL,OUTCONF,OPTIONSCELL) performs the Unscented   
%  Kalman Filter algorithm on a Discrete State-Space object using  
%  data from the simulation model SIMMODEL.
%  The parameters OUTCONF and OPTIONSCELL can be omitted if desired but
%  the order of the parameters must remain.
%  For example: UKF(SSOBJ,SIMMODEL,OPTIONSCELL) and UKF(SSOBJ,Y,OUTCONF) are
%  allowed, but UKF(SSOBJ,SIMMODEL,OPTIONSCELL,OUTCONF) is not allowed
%
% Reference: Optimal State Estimation - Chapter 14 - Dan Simon 

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
            error('DA:StateSpaceModels:ss_D:UKF:cellDim','Cell must be vector sized.')
        else
            ni=1+max(size(varargin,1),size(varargin,2));
        end
	end    
    
    % Retrieve k0 from ssObj       
    k0=ssObj.k0;
    
    % Unscented Kalman Filter should be applied on additive noise models
    if isa(ssObj,'ss_DNL') || isa(ssObj,'ss_DNL_AMN')
        warning('DA:StateSpaceModels:ss_D:UKF:ssWarning','Unscented Kalman Filter should only be applied on additive noise models.')
    end    
    
    % Unscented Kalman Filter should be applied on zero mean gaussian noise       
    if ~iszeromean(mean(ssObj.w,k0)) || ~iszeromean(mean(ssObj.v,k0)) || ~isa(ssObj.w,'nm_gauss') || ~isa(ssObj.v,'nm_gauss')
        warning('DA:StateSpaceModels:ss_D:UKF:nonZeroMean','Noise models should be zero mean and Gaussian for UKF.')
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
        error('DA:StateSpaceModels:ss_D:UKF:classMismatch','Class mismatch for y: expected sim_D object or double.')
    end

    %Resolve optionsCell
    if ~isempty(optionsCell)
        ni=max([size(optionsCell,1) size(optionsCell,2)]);
        if ni>0
            n=size(optionsCell{1});
            if isa(optionsCell{1},'double')&&length(n)==2&&n(1)==1&&n(2)==1         %scalar
                kappa=optionsCell{1};    
            	if isempty(kappa);kappa=0;end
            else
                error('DA:StateSpaceModels:ss_D:UKF:kappa','Kappa must be a scalar.')
            end
        else
            kappa=0;             
        end
        if ni>1
            n=size(optionsCell{2});
            if isa(optionsCell{2},'double')&&length(n)==2&&n(1)==1&&n(2)==1         %scalar
                spTwice=optionsCell{2};    
            	if isempty(spTwice);spTwice=0;end  
            else
                error('DA:StateSpaceModels:ss_D:UKF:spTwice','SpTwice must be a scalar.')
            end
        else
            spTwice=0;  
        end
    else
        kappa=0;
        spTwice=0;        
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
        warning('DA:StateSpaceModels:ss_D:UKF:uSize','Amount of inputs u is not consistent with amount of measurements y. Smallest amount is used.')
        if size(u,2)<size(y,2);samples=size(u,2);end
        if size(u,2)>size(y,2);samples=size(y,2);end
    end        
    
    %Timings        
	Ts=ssObj.Ts;
	kIndex=k0:1:(k0+samples-1);	

    %========================================%
    %   UNSCENTED KALMAN FILTER ALGORITHM    %
    %========================================%
    
    %Initial estimates
    xhat = mean(ssObj.x0,k0);
    Phat = cov(ssObj.x0,k0);
    nStates=size(xhat,1);
	nOutputs=size(y,1);
    %pre-allocate memory of data arrays & fill in first sample
    if xfRet; xfDa = zeros(nStates,samples);xfDa(:,1)=xhat;end;
    if ~xfRet;xfDa=[];end;
    xaDa = zeros(nStates,samples);xaDa(:,1)=xhat;
    if PfRet && covFull; PfDa = zeros(nStates,nStates,samples);PfDa(:,:,1)=Phat; end;
    if PfRet && ~covFull; PfDa = zeros(nStates,samples);PfDa(:,1)=diag(Phat); end;
    if ~PfRet;PfDa=[];end;
    if PaRet && covFull; PaDa = zeros(nStates,nStates,samples);PaDa(:,:,1)=Phat;end;
    if PaRet && ~covFull; PaDa = zeros(nStates,samples);PaDa(:,1)=diag(Phat);end;
    if ~PaRet;PaDa=[];end;
    
    if (nStates+kappa)==0
        error('DA:StateSpaceModels:ss_D:UKF:kappa','Amount of states and kappa cannot sum to 0.')
    end
        
    W0=kappa/(nStates+kappa);
    Wi=0.5/(nStates+kappa) + zeros(1,2*nStates);
    W=[W0 Wi];
    
    for i=2:samples

        %FORECAST STEP
        xSp=sigmaPoints(xhat,Phat,nStates+kappa); 
        [xhat,xSp,Phat,xDev]=UT(ssObj,'f',xSp,W,nStates,kIndex(i-1),getU(u,i-1));
        if xfRet; xfDa(:,i)=xhat;end;                       %SAVE xf
        if PfRet && covFull; PfDa(:,:,i)=Phat;end;          %SAVE Pf
        if PfRet && ~covFull; PfDa(:,i)=diag(Phat);end;     %SAVE Pf        
        
        %MEASUREMENT UPDATE
        if spTwice
            xSp=sigmaPoints(xhat,Phat,nStates+kappa); 
        end
        [yhat,~,Pyy,yDev]=UT(ssObj,'h',xSp,W,nOutputs,kIndex(i),getU(u,i));
        Pxy=xDev*diag(W)*yDev';
        
        %KALMAN GAIN
        K = Pxy / (Pyy);          
        
        %ANALYSIS STEP
        xhat = xhat + K*( y(:,i)-yhat );    %xa
        Phat = Phat - K * Pyy * K';         %Pa
        Phat = (Phat + Phat')/2; %enforce symmetry           
        xaDa(:,i)=xhat;                                     %SAVE xa        
        if PaRet && covFull; PaDa(:,:,i)=Phat;end;          %SAVE Pa
        if PaRet && ~covFull; PaDa(:,i)=diag(Phat);end;  	%SAVE Pa           
        
    end
      
    %MAKE UNREQUIRED DATA EMPTY    
	if ~uRet;u=[];end;
    if ~yRet;y=[];end;
    if ~xtrueRet;xtrue=[];end;
        
    if ~isempty(u);u=u(:,1:samples);end;
    if ~isempty(y);y=y(:,1:samples);end;
    if xtrueRet && simTrue;xtrue=simObj.x(:,1:samples);end;    

    % + to change to double
    daModel = dam_D('UKF',xtrue,xfDa,xaDa,PfDa,PaDa,u ,y ,kIndex,+(covFull),Ts,ssObj.TimeUnit);  
                    
end
            
function xSp=sigmaPoints(x,P,nPlusKappa)
% Sigma points around reference point
% INPUTS
% x:             state vector
% P:             covariance matrix
% nPlusKappa: 	 nStates+kappa
%
% OUTPUTS
% xSp:           sigma points states

    xTi = sqrt(nPlusKappa)*chol(P)';
    xi = repmat(x,1,length(x));
    xSp = [x xi+xTi xi-xTi]; 

end

function [trMean,trSp,trCov,dev]=UT(ssObj,func,sp,W,nOut,t,u)
% Unscented Transformation
% INPUTS
% func:      function f or h
% sp:        sigma points
% W:         weights
% nOut:      output size
% t:         current time
% u:         inputs
%
% OUTPUTS
% trMean:    transformed mean
% trSp:      transformed sigma points
% trCov:     transformed covariance
% dev:       deviation transformed sigmapoints - transformed mean

    nSp=size(sp,2);
    trMean=zeros(nOut,1);
    trSp=zeros(nOut,nSp);
        
    if strcmp(func,'f')
        covN = cov(ssObj.w,t);
    elseif strcmp(func,'h')
        covN = cov(ssObj.v,t);      
    end
    
    for i=1:nSp 
        if strcmp(func,'f')
          	muW=mean(ssObj.w,t);   
            trSp(:,i)=eval_ftot(ssObj,sp(:,i),t,u,muW); 
        elseif strcmp(func,'h')
          	muV=mean(ssObj.v,t);   
            trSp(:,i)=eval_htot(ssObj,sp(:,i),t,u,muV);             
        end
        trMean=trMean+W(i)*trSp(:,i);       
    end
    
    dev=trSp-trMean(:,ones(1,nSp));
    trCov=dev*diag(W)*dev'+covN;   
    
end
    
    
    
    
  