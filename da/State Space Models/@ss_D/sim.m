function simModel=sim(ssObj,samples,u,conf,x0,noise)
%SIM simulates any discrete State-Space object
%
%  - Input variable(s) -
%  SSOBJ: any discrete State-Space object. (type 'help ss_D')
%
%  SAMPLES: amount of samples/iterations to be simulated. Default=1.
%
%  U: matrix of inputs where every column is a input vector. Can be left 
%  empty ([ ]) if desired. 
%
%  CONF: string that indicates which variables must reside in the
%  simulation model. Possible string values contain 'x', 'u' and 'y'. 
%  For example: 'x y' retains x and y matrices.
%  Default: 'x u y'
%
%  X0: initial state column vector. If left empty ([ ]), a sample will be 
%  taken from the x0 noise model of the State-Space model.  
%
%  NOISE: option to specify noise behaviour during simulation.
%   * 0: no noise
%   * 1: only simulate process noise
%   * 2: only simulate measurement noise
%   * 3: simulate process and measurement noise (Default)
%
%  - Output variable(s) -
%  SIMMODEL: a discrete time State-Space simulation model of type sim_D.
%  
%  - Construction -          
%  SIMMODEL=SIM(SSOBJ,SAMPLES,U,CONF,X0,NOISE) simulates the discrete 
%  State-Space object SSOBJ.
%
%  SIMMODEL=SIM(SSOBJ,SAMPLES,U,CONF,X0) simulates the discrete State-Space 
%  object SSOBJ with both process and measurement noise.
%
%  SIMMODEL=SIM(SSOBJ,SAMPLES,U,CONF) simulates the discrete State-Space 
%  object SSOBJ with both process and measurement noise and with the 
%  initial state X0 a sample from the x0 noise model of SSOBJ.
%
%  SIMMODEL=SIM(SSOBJ,SAMPLES,U) simulates the discrete State-Space 
%  object SSOBJ with both process and measurement noise and with the 
%  initial state X0 a sample from the x0 noise model of SSOBJ, CONF='x u y'.
%
%  SIMMODEL=SIM(SSOBJ,SAMPLES) simulates the discrete State-Space 
%  object SSOBJ with both process and measurement noise and with the 
%  initial state X0 a sample from the x0 noise model of SSOBJ, CONF='x u y'
%  and no inputs (U=[]).
%
%  SIMMODEL=SIM(SSOBJ) simulates one sample of the discrete State-Space 
%  object SSOBJ with both process and measurement noise and with the 
%  initial state X0 a sample from the x0 noise model of SSOBJ, CONF='x u y'
%  and no inputs (U=[]).
%

    % check amount of input arguments
	narginchk(1,6);
	ni = nargin;  

	if ni<2;samples=1;end;
	if ni<3;u=[];end;
	if ni<4;conf='x u y';end;
	if ni<5;x0=[];end;    
    if ni<6;noise=3;end;
            
	if isempty(samples);samples=1;end;
	if isempty(conf);conf='x u y';end;
    if isempty(noise);noise=3;end;
            
    % check content of input arguments
	samples=simCheckArgs(samples,'samples');
	u=simCheckArgs(u,'u');
	conf=simCheckArgs(conf,'conf');
	x0=simCheckArgs(x0,'x0');    
	noise=simCheckArgs(noise,'noise');        
           
    % determine which variables must be retained
	xRet = ~isempty(strfind(conf,'x'));
	uRet = ~isempty(strfind(conf,'u'));
            
	%samples can not be larger than amount of u's
	if ~isempty(u) && (samples > size(u,2))
        samples=size(u,2);
	end
    %obtain k0 and Ts from ssobj        
	k0=ssObj.k0;
	Ts=ssObj.Ts;
	kIndex=k0:1:(k0+samples-1);	
            
	%acquire x0 from input of from noise model within ssobj
    if isempty(x0)
        x = sample(ssObj.x0,1,k0); 
    else
        x=x0;
    end
    
    %noise int to bits
    switch noise
        case 0
        pn=0;mn=0;    
        case 1
        pn=1;mn=0;            
        case 2
        pn=0;mn=1;            
        case 3
        pn=1;mn=1;            
    end
    
    %already evaluate first vars to determine sizes
	if mn==1;y = eval_htot(ssObj,x,k0,getU(u,1),[]);end   %y0
	if mn==0;y = eval_htot(ssObj,x,k0,getU(u,1),0);end   %y0    
       
	%pre-allocate memory & fill in first sample for x and y
	if xRet; xSim = zeros(size(x,1),samples);xSim(:,1)=x;end;
	ySim = zeros(size(y,1),samples);ySim(:,1)=y;          

    %obtain other samples
	for i=2:samples
        %states       
        if pn==1;x = eval_ftot(ssObj,x,kIndex(i-1),getU(u,i-1),[]);end
        if pn==0;x = eval_ftot(ssObj,x,kIndex(i-1),getU(u,i-1),0);end
        %measurements
        if mn==1;y = eval_htot(ssObj,x,kIndex(i),getU(u,i),[]);end;
        if mn==0;y = eval_htot(ssObj,x,kIndex(i),getU(u,i),0);end;
        
        if xRet; xSim(:,i)=x;end;
        ySim(:,i)=y; 
                
	end
            
    %make simModel
    if ~uRet;u=[];end;
    if ~xRet;xSim=[];end;
    
    if ~isempty(u);u=u(:,1:samples);end;
    
    simModel = sim_D(xSim,u,ySim,kIndex,Ts,ssObj.TimeUnit);
            
end