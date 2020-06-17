% Demo about creating and manimulating a nonlinear state space model
% ss_DNL

clc;

disp('This is a demo of how to create and manipulate a nonlinear ');
disp('state space model (ss_DNL). It is assumed that the demo');
disp('''demo_ss_DL'' has already been studied by the user. ');
disp('Note that this demo is very similar to the demo ''demo_ss_DNL_AN''.');
disp('It has been added for completeness sake and to illustrate the similarities.');
disp('The model equations of a ss_DNL model are:');
disp(' ');
disp(' x(k+1) = f(x(k),u(k),w(k),k)');
disp(' y(k)   = h(x(k),u(k),v(k),k)');
disp(' ');
disp('So, to configure a  model the following');
disp('parameters are required:');
disp('- The function handles for f and h');
disp('- The process and measurement noise models w and v');
disp('- The intial state estimate x0');
disp('where x0, w and v can be any kind of noise model and the system matrices');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('We start with a simple system of two states which are both measured');
disp('and and there is one input. ');
disp('The function handles for f is @demo_DNL_f and is configured as');
disp(' ');
disp('function x1 = demo_DNL_f(x,t,u,w,Ts)');
disp(' ');
disp('	x1(1,1) = x(1)+2*x(2) + u + w(1);');
disp('	x1(2,1) = x(1)+x(2) + w(2);');
disp(' ');
disp('end');
disp(' ');
disp('The function handle for h is @demo_DNL_h and is configured as');
disp(' ');
disp('function y = demo_DNL_h(x,t,u,v,Ts)');
disp(' ');
disp('	y(1,1) = x(1) + v(1);');
disp('	y(2,1) = x(2) + v(2);');
disp(' ');
disp('end');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('For simplicity the noise models will be configured as LTI noise');
disp('models, but LTV and function handle are of course also possible.');
disp('Let us define the process and measurement noise as zero mean');
disp('uncorrelated Gaussian noise and the initial state x0 is defined ');
disp('with an initial state [1;2] and unit covariance');
disp(' ');
disp('>> w = nm_gauss_lti([0;0],[2 0;0 2]);');
disp('>> v = nm_gauss_lti([1 0;0 1]);');
disp('>> x0 = nm_gauss_lti([1;2];');
w = nm_gauss_lti([0;0],[2 0;0 2]);
v = nm_gauss_lti([1 0;0 1]);
x0 = nm_gauss_lti([1;2]);
disp(' ');
disp('Now that the function handles and the noise sources are defined we are');
disp('ready to configure the state space model.');
disp(' ');
disp('>>ssObj=ss_DNL(@demo_DNL_f,@demo_DNL_h,x0,w,v);');
ssObj=ss_DNL(@demo_DNL_f,@demo_DNL_h,x0,w,v);
disp(' ');
disp('The more general way to generate a DNL state space model is ');
disp('briefly explained now to summarize all possible arguments. The most');
disp('general form is:');
disp(' ');
disp(' ''ss_DNL(@demo_DNL_f,@demo_DNL_h,x0,w,v,k0,Ts,TimeUnit,...');
disp('          ,@demo_DNL_fJacX,@demo_DNL_fJacW,@demo_DNL_hJacX,@demo_DNL_hJacV ) ''');
disp(' ');
disp('where k0 is the initial step number, Ts the sample time and TimeUnit');
disp('defines the unit of the time. The last four parameters are function');
disp('handles that return the Jacobians of f with respect to x and w, and ');
disp('of h with respect to x and v. We will come back to this later.');
disp('For now, the obtained state space object is sufficient.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Lets see how this state space model is presented.');
disp(' ');
disp('>>ssObj');
ssObj
disp('The parameters k0, Ts and TimeUnit contain their default values and');
disp('the Jacobians are undefined, as expected.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now, the different methods of this class are summarized.');
disp('A first method allows to check if the model is properly configured,');
disp('i.e. that the functions are properly returning values. This is done');
disp('with the function ''checkConsist'' which requires the input size and');
disp('the measurement size as extra input arguments and returns 1 if all');
disp('which is of course the case in this demo.');
disp(' ');
disp('>>uSize=1; ySize=2; ');
disp('>>checkConsist(ssObj,uSize,ySize)');
uSize=1; ySize=2;
checkConsist(ssObj,uSize,ySize)
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('A second method enables to retrieve a state or measurement update,');
disp('where the state equation is ftot and the measurement equation htot.');
disp('In this case ftot and htot are the same as f and h.');
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
disp('This short-hand notation can not be applied to u.');
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
disp('system equations ftot and htot.');
disp(' ');
disp('There are four possible Jacobian functions: ');
disp('- eval_ftotJacX : Calculates the Jacobian of ftot with respect to x.');
disp('                  Note that this is the same as f with respect to x.');
disp('- eval_htotJacX : Calculates the Jacobian of htot with respect to x.');
disp('                  Note that this is the same as h with respect to x.');
disp('- eval_ftotJacW : Calculates the Jacobian of ftot with respect to w.');
disp('                  Note that this is the same as f with respect to w.');
disp('- eval_htotJacV : Calculates the Jacobian of htot with respect to v.');
disp('                  Note that this is the same as h with respect to v.');
disp(' ');
disp('All of them have an identical syntax, which is identical to the one');
disp('of eval_ftot. Here, the method eval_ftotJacX is demonstrated: ');
disp(' ');
disp('>>eval_ftotJacX(ssObj,x1,k1,u1,w1)');
eval_ftotJacX(ssObj,x1,k1,u1,w1)
disp(' ');
disp('This is correct, even without supplying any Jacobian matrices. When');
disp('no Jacobian equation is provided, the toolbox will estimate the Jacobian');
disp('numerically.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('To provide the Jacobian to the toolbox, the following function is used');
disp(' ');
disp('function JacX = demo_DNL_fJacX(x,k,u,~,Ts)');
disp(' ');
disp('	JacX = [1 2 ; 1 1] ;');
disp(' ');
disp('end');
disp(' ');
disp('and it is provided to the model as follows');
disp(' ');
disp('>>ssObj.fJacX=@demo_DNL_fJacX;');
ssObj.fJacX=@demo_DNL_fJacX;
disp(' ');
disp('now the result is');
disp(' ');
disp('>>eval_ftotJacX(ssObj,x1,k1,u1,w1)');
eval_ftotJacX(ssObj,x1,k1,u1,w1)
disp(' ');
disp('which is of course again correct.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('And to end with some more Jacobians ...');
disp(' ');
disp('>>eval_ftotJacW(ssObj,x1,k1,u1,w1)');
eval_ftotJacW(ssObj,x1,k1,u1,w1)
disp(' ');
disp('>>eval_htotJacX(ssObj,x1,k1,u1,w1)');
eval_htotJacX(ssObj,x1,k1,u1,w1)
disp(' ');
disp('This concludes this demo.');

