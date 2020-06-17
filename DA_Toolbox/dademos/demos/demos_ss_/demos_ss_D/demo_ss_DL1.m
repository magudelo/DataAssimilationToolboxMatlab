% Demo about creating and manimulating a linear time variant state space model
% ss_DL

clc;

disp('This is a demo of how to create and manipulate a linear time variant ');
disp('state space model(ss_DL).');
disp('Note that in its most general form the model equations are:');
disp(' ');
disp(' x(k+1)  = A(k)x(k) + B(k)u(k) + w(k) ');
disp(' y(k)    = C(k)x(k) + D(k)u(k) + v(k) ');
disp(' ');
disp('So, to configure a linear time variant state space model the following');
disp('parameters are required:');
disp('- The system matrices A, B, C and D');
disp('- The process and measurement noise models w and v');
disp('- The intial state estimate x0');
disp('where x0, w and v can be any kind of noise model and the system matrices');
disp('A, B, C and D which can be either LTI (matrix) or LTV (3D-array).');
disp(' ');
disp('Hence, a large amount of configurations are possible. Let us start with');
disp('LTI system matrices and focus on how to enter the noise models. ');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('We start with a simple system of two states which are both measured');
disp('and and there is one input. The system matrices are defined as');
disp(' ');
disp('>> A = [1 2;1 1];');
disp('>> B = [1 ; 0];');
disp('>> C = [1 0;0 1];');
disp('>> D = 0;');
disp(' ');
A = [1 2;1 1];
B = [1 ; 0];
C = [1 0;0 1];
D = 0;
disp('Note that since there is no feedthrough, the matrix D can be set');
disp('to zero, the toolbox will properly resize it.');
disp(' ');
disp('The noise models can be entered in several intuitive ways:');
disp('- either first configure a noise model and use it as input argument');
disp('- enter the noise model as a cell array');
disp('- enter either the mean mu or the covariance Sigma');
disp('When using the last two methods, the noise is considered to be');
disp('Gaussian. The toolbox will automatically decide whether the noise');
disp('model is LTI, LTV or of type function handle. See the help files or ');
disp('demos of the Gaussian noise models for the different techniques to');
disp('configure such a noise model.');
disp('For simplicity the noise models will be configured as LTI noise');
disp('models, but LTV and function handle are of course also possible.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Let us already define the process and measurement noise as zero');
disp('mean uncorrelated Gaussian noise. Note that the different configuration');
disp('styles both give zero mean noise.');
disp(' ');
disp('>> w = nm_gauss_lti([0;0],[2 0;0 2]);');
disp('>> w');
w = nm_gauss_lti([0;0],[2 0;0 2]);
w
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Let us already define the process and measurement noise as zero');
disp('mean uncorrelated Gaussian noise. Note that the different configuration');
disp('styles both give zero mean noise.');
disp(' ');
disp('>> v = nm_gauss_lti([1 0;0 1]);');
disp('>> v');
v = nm_gauss_lti([1 0;0 1]);
v
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that the system matrices and the noise sources are defined we are');
disp('ready to configure the state space model in three different ways');
disp('that provide the same result. First, lets define x0 explicitly');
disp('with initial state [1;2] and unit covariance');
disp(' ');
disp('>> x0 = nm_gauss_lti([1;2],[1 0;0 1]);');
x0 = nm_gauss_lti([1;2],[1 0;0 1]);
disp(' ');
disp('and construct the ss_DL model as');
disp(' ');
disp('>>ssObj=ss_DL(A,B,C,D,x0,w,v);');
ssObj=ss_DL(A,B,C,D,x0,w,v);
disp(' ');
disp('To see how the noise model is configured you can type');
disp(' ');
disp('>>ssObj.x0');
ssObj.x0
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Lets try the second way, using a cell array, to configure the noise');
disp('model in a ss_D object');
disp(' ');
disp('>>ssObj=ss_DL(A,B,C,D,{[1;2],[1 0;0 1]},w,v);');
ssObj=ss_DL(A,B,C,D,{[1;2],[1 0;0 1]},w,v);
disp(' ');
disp('and check if it is configured properly');
disp(' ');
disp('>>ssObj.x0');
ssObj.x0
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Since the unit covariance is the default Sigma for a Gaussian noise');
disp('model, we do not have to define it explicitly, the mean is sufficient');
disp(' ');
disp('>>ssObj=ss_DL(A,B,C,D,[1;2],w,v);');
ssObj=ss_DL(A,B,C,D,[1;2],w,v);
disp(' ');
disp('and check if it is configured properly');
disp(' ');
disp('>>ssObj.x0');
ssObj.x0
disp(' ');
disp('Note that, although allowed, the cell array was not used for one input');
disp('argument. Also notice that these three ways of entering a noise model');
disp('can also be applied on the w and v.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('To obtain an LTV state space model just define the LTV system matrices');
disp('as 3D arrays. The priciple of obtaining the correct matrix from a ');
disp('3D array is identical as for a LTV Gaussian noise model. As such, ');
disp('this is not included in this demo. Instead run through the demo ');
disp('''demo_nm_gauss_ltv1'' or consult the help file of ''ss_DL''');
disp('for more information.');
disp(' ');
disp('Still, the most general way to generate a state space model is ');
disp('briefly explained now to summarize all possible arguments. The most');
disp('general form is:');
disp(' ');
disp(' ''ss_DL(A,B,C,D,x0,w,v,k0,Ts,TimeUnit,kIndex,kMethod) ''');
disp(' ');
disp('where k0 is the initial step number, Ts the sample time, TimeUnit');
disp('defines the unit of the time and kIndex and kMathod define how');
disp('information is extracted from the system matrices.');
disp('So, lets change some settings as follows:');
disp(' ');
disp('>>ssObj=ss_DL(A,B,C,D,x0,w,v,2,10,''minutes'',[],''high'');');
ssObj=ss_DL(A,B,C,D,x0,w,v,2,10,'minutes',[],'high');
disp(' ');
disp('which makes the initial step number two and the sampling time is now');
disp('10 minutes. The last two arguments have no meaning in this demo since');
disp('none of the system matrices are time variant.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Lets see how the state space model is presented.');
disp(' ');
disp('>>ssObj');
ssObj
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('To change any property just proceed as follows:');
disp(' ');
disp('>>ssObj.k0=0');
ssObj.k0=0
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('To conclude this demo the different methods are summarized.');
disp('A first method enables to retrieve the system matrices at a');
disp('certain step number. For instance the matrix A at step k=5: ');
disp(' ');
disp('>>k=5;');
disp('>>eval_A(ssObj,k)');
k=5;
eval_A(ssObj,5)
disp(' ');
disp('Of course, since matrix A is LTI, the step number does not need to ');
disp('be supplied. Still, it it recommended to always use the most general');
disp('notation so that your code does not need alterations when a system ');
disp('matrix is changed from LTI to LTV. Do not call the property A for');
disp('the same reas! The other matrices are called with eval_B, eval_C');
disp('and eval_D.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('A second method enables to retrieve a state or measurement update,');
disp('where the state equation is ftot and the measurement equation htot.');
disp('Do not confuse this with f and h, ftot/htot are the entire equations!');
disp('To obtain a state update from current state x1, at step k1, with input u1');
disp('and w1=[0;0] the following syntax is used (for all kinds of ss_D models)');
disp(' ');
disp('>>k1=1; x1=[1;1]; u1=2; w1=[0;0]; ');
disp('>>eval_ftot(ssObj,x1,k1,u1,w1)');
k1=1;x1=[1;1]; u1=2; w1=[0;0];
eval_ftot(ssObj,x1,k1,u1,w1)
disp(' ');
disp('the notation can be made shorter for w. If the zero vector is required');
disp('it is sufficient to use 0 which provides the same result:');
disp(' ');
disp('>>eval_ftot(ssObj,x1,k1,u1,0)');
eval_ftot(ssObj,x1,k1,u1,0)
disp(' ');
disp('The same is valid for u.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Another common procedure is to have a state update using the expected');
disp('value of the noise source w as follows');
disp(' ');
disp('>>eval_ftot(ssObj,x1,k1,u1,mean(ssObj.w,k1))');
eval_ftot(ssObj,x1,k1,u1,mean(ssObj.w,k1))
disp(' ');
disp('which gives the same result since w is zero mean in this case.');
disp('It is very important to remember that a sample is taken from the ');
disp('noise source by using an empty matrix');
disp(' ');
disp('>>eval_ftot(ssObj,x1,k1,u1,[])');
eval_ftot(ssObj,x1,k1,u1,[])
disp(' ');
disp('The method to update the measurement equation is analogous with');
disp('syntax ''eval_htot(ssObj,x1,k1,u1,v1)''.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('The last range of methods are used to obtain the Jacobians of the ');
disp('equations. Of course, in case of this linear model, the ');
disp('Jacobians are equal to their corresponding system matrices. Still,');
disp('they are implemented to allow a generic interface to both linear');
disp('and nonlinear state space models.');
disp(' ');
disp('There are four possible Jacobian functions: ');
disp('- eval_ftotJacX : Calculates the Jacobian of ftot with respect to x.');
disp('                  This is the matrix A for this model.');
disp('- eval_htotJacX : Calculates the Jacobian of htot with respect to x.');
disp('                  This is the matrix C for this model.');
disp('- eval_ftotJacW : Calculates the Jacobian of ftot with respect to w.');
disp('                  This is a unity matrix for this model.');
disp('- eval_htotJacV : Calculates the Jacobian of htot with respect to v.');
disp('                  This is a unity matrix for this model.');
disp(' ');
disp('All of them have an identical syntax, which is identical to the one');
disp('of eval_ftot. Here, the method eval_ftotJacX is demonstrated: ');
disp(' ');
disp('>>eval_ftotJacX(ssObj,x1,k1,u1,w1)');
eval_ftotJacX(ssObj,x1,k1,u1,w1)
disp(' ');
disp('Of course, as mentioned before, this yields the system matrix A.');
disp(' ');
disp('This concludes this demo.');

