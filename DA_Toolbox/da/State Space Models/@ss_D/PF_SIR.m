function daModel=PF_SIR(ssObj,varargin)
%PF_SIR performs the Sampling Importance Resampling (SIR) Filter algorithm  
%    on a Discrete State-Space object
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
%  data assimilation model. Possible string values contain 'Pa', 'xF', 'Pf' 
%  ,'y', 'u', 'cov', 'var'. 'xa' is always returned.
%  For example: 'Pa y' retains xa, Pa and y matrices.
%  Default: 'xa'
%
%  OPTIONSCELL: cell that contains additional options. 
%  The first variable has the option NPARTICLES or INITPARTICLES. 
%  NPARTICLES indicates the number of particles to use. When NPARTICLES is 
%  used, the particles are selected randomly from the X0 noise model of the 
%  SSOBJ. INITPARTICLES contains user defined initial particles. When 
%  neither NPARTICLES or INITPARTICLES are configured NPARTICLES is set 
%  equal to 100.
%  The second variable is the string RESAMPLER which indicates the resample
%  algorithm. Possible values are:
%   * 'multinom':       see reference [3]
%   * 'residual':       see reference [3]
%   * 'systematic':     see reference [3]
%   * 'stratisfied':    see reference [3] (Default)
%  When the second variable is selected as 'residual' a third variable can
%  be added: 'multinom', 'systematic' or 'stratisfied'.
%
%  - Output variable(s) -
%  DAMODEL: a discrete time State-Space data assimilation model of type dam_D.
%  
%  - Construction -          
%  DAMODEL=PF_SIR(SSOBJ,Y,U,OUTCONF,OPTIONSCELL) performs the SIR Particle
%  Filter algorithm on a Discrete State-Space object.
%  The parameters U, OUTCONF and OPTIONSCELL can be omitted if desired but
%  the order of the parameters must remain.
%  For example: PF_SIR(SSOBJ,Y,U,OPTIONSCELL) and PF_SIR(SSOBJ,Y,OUTCONF) are
%  allowed, but PF_SIR(SSOBJ,OUTCONF,Y) is not allowed
%
%  DAMODEL=PF_SIR(SSOBJ,SIMMODEL,OUTCONF,OPTIONSCELL) performs the SIR  
%  Particle Filter algorithm on a Discrete State-Space object using  
%  data from the simulation model SIMMODEL.
%  The parameters OUTCONF and OPTIONSCELL can be omitted if desired but
%  the order of the parameters must remain.
%  For example: PF_SIR(SSOBJ,SIMMODEL,OPTIONSCELL) and PF_SIR(SSOBJ,Y,OUTCONF) 
%  are allowed, but PF_SIR(SSOBJ,SIMMODEL,OPTIONSCELL,OUTCONF) is not allowed
%
% References: 
% [1] A Tutorial on Particle Filters for Online Nonlinear/Non-Gaussian 
% Bayesian Tracking Arulampalam et al - IEEE TRANSACTIONS ON SIGNAL 
% PROCESSING, VOL. 50, NO. 2, FEBRUARY 2002
%
% [2] Resampling in particle filters - Internship performed at Division of 
% Automatic Control Department of Electrical Engineering Linkopings 
% University, Sweden - Jeroen D. Hol

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
            error('DA:StateSpaceModels:ss_D:PF_SIR:cellDim','Cell must be vector sized.')
        else
            ni=1+max(size(varargin,1),size(varargin,2));
        end
	end    
    
    % Retrieve k0 from ssObj       
    k0=ssObj.k0;
    
    % SIR Particle Filter should be applied on additive measurement noise models
    if isa(ssObj,'ss_DNL')
        warning('DA:StateSpaceModels:ss_D:PF_SIR:ssWarning','SIR Particle Filter should only be applied on additive measurement noise models.')
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
        error('DA:StateSpaceModels:ss_D:PF_SIR:classMismatch','Class mismatch for y: expected sim_D object or double.')
    end

    %Resolve optionsCell
    if ~isempty(optionsCell)
        ni=max([size(optionsCell,1) size(optionsCell,2)]);
        if ni>0
            n=size(optionsCell{1});
            if isa(optionsCell{1},'double')&&length(n)==2&&n(1)==1&&n(2)==1         %scalar
                nParticles=optionsCell{1};
                initParticles=0;
            elseif isa(optionsCell{1},'double')&&length(n)==2&&n(1)>1&&n(2)>1       %matrix 
                initParticles=optionsCell{1};
                nParticles=size(initParticles,2);
                if size(initParticles,1)~=ssObj.x0
                    error('DA:StateSpaceModels:ss_D:PF_SIR:initParticles','Number of states of initParticles is incorrect.')
                end  
            elseif isempty(optionsCell{1})
                nParticles = 100; 
                initParticles=0;
            else
                error('DA:StateSpaceModels:ss_D:PF_SIR:argMismatch','Optionscell must contain scalar nParticles or initParticles as first argument.')
            end  
        else
            nParticles = 100;  
            initParticles=0;
        end   
        if ni>1
            if isa(optionsCell{2},'char')                                           %string
                resampler=optionsCell{2};
            elseif isempty(optionsCell{2})
                resampler = 'stratisfied';                
            else
                error('DA:StateSpaceModels:ss_D:PF_SIR:argMismatch','Optionscell must contain Resampler string as second argument.')
            end     
        else
            resampler = 'stratisfied'; 
        end         
        if strcmp(resampler,'residual')&&ni>2
            if isa(optionsCell{3},'char')                                           %string
                resampler2=optionsCell{3};
            elseif isempty(optionsCell{3})
                resampler2 = 'stratisfied';                
            else
                error('DA:StateSpaceModels:ss_D:PF_SIR:argMismatch','Optionscell must contain Resampler string as third argument.')
            end     
        else
            resampler2 = 'stratisfied'; 
        end                 
    else
        nParticles = 100;   
        initParticles=0;
        resampler = 'stratisfied'; 
        resampler2 = 'stratisfied'; 
    end 
    if nParticles<2
        error('DA:StateSpaceModels:ss_D:PF_SIR:nParticles','Number of particles must be minimal equal to two.')
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
        warning('DA:StateSpaceModels:ss_D:PF_SIR:uSize','Amount of inputs u is not consistent with amount of measurements y. Smallest amount is used.')
        if size(u,2)<size(y,2);samples=size(u,2);end
        if size(u,2)>size(y,2);samples=size(y,2);end
    end        
    
    %Timings        
	Ts=ssObj.Ts;
	kIndex=k0:1:(k0+samples-1);	

    %=========================================%
    %      SIR PARTICLE FILTER ALGORITHM      %
    %=========================================%
    
    %Initial particles: if empty create random particles    
    if initParticles==0 || isempty(initParticles)
        initParticles=sample(ssObj.x0,nParticles,k0);                  %initial particles
    end        

    %Initial estimates & pre-allocate memory
    xhat_par=initParticles;
    xhat = mean(xhat_par,2);
    nStates=size(xhat,1);
	nOutputs=size(y,1);
    yf_par=zeros(nOutputs,nParticles);
    resampled=1;                                                    %initial behaviour is identical to resampled
    if PfRet || PaRet
        Ex=xhat_par-repmat(xhat,1,nParticles);                    	%particle error covariance matrix
        [Phat] = P_Save_Stats(Ex,covFull);
    end
    wi=ones(nParticles,1) / nParticles;                             %initial weights column vector
    
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
        for j=1:nParticles  %forecast every particle  
            xhat_par(:,j) = eval_ftot(ssObj,xhat_par(:,j),kIndex(i-1),getU(u,i-1),[]);   % xf_par(nxN)
            yf_par(:,j) = eval_htot(ssObj,xhat_par(:,j),kIndex(i),getU(u,i),muV);        % (mxN)
        end  
        
        %SAVE STATS
        [xhat,Phat] = PF_Save_Stats(xhat_par,wi,xfRet,PfRet,covFull,resampled);        
        if xfRet; xfDa(:,i)=xhat;end;                     %SAVE xf
        if PfRet && covFull; PfDa(:,:,i)=Phat;end;        %SAVE Pf
        if PfRet && ~covFull; PfDa(:,i)=Phat;end;         %SAVE Pf
        
        %ANALYSIS STEP        
        %Assign the particle the weights w_i
        wi=pdf(ssObj.v,repmat(y(:,i),1,nParticles)-yf_par,kIndex(i));	% Calculate weights
        sumwi=sum(wi);
        if sumwi <= realmin
            error('DA:StateSpaceModels:ss_D:PF_SIR:sumWi','The sum of weights is zero. Please adjust the parameters or the method.')
        end
        wi=wi/sumwi;                                                        % Normalize weights
	
        %resample
        resamplerHandle=str2func(strcat('PFRS_',resampler));
        [xhat_par,wi]=feval(resamplerHandle,xhat_par,wi,length(wi),resampler2);  % Resample 
        resampled=1;      

        %SAVE STATS              
        [xhat,Phat] = PF_Save_Stats(xhat_par,wi,1,PaRet,covFull,resampled);       
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
    daModel = dam_D('PF_SIR',xtrue,xfDa,xaDa,PfDa,PaDa,u ,y ,kIndex,+(covFull),Ts,ssObj.TimeUnit);  
                         
end

    
    
  