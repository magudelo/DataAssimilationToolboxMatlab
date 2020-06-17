% Demo about using the simulation method to create simulation objects.

clc;

disp('This is a demo of how to apply the simulation model on a discrete-time ');
disp('state space model to obtain a simulation model.');
disp('This demo uses the preconfigured functions ''demo_sim1f'' and ''demo_sim1h''');
disp('which represent the Van Der Pol oscillator.');
disp(' ');
disp('First a ss_DNL_AN state space model is created using the following steps:');
disp(' ');
disp('1) Create the noise models ');
disp('>> mu_x0=[0 5]'';sigma_x0=diag([5 5]);');
disp('>> mu_w=[0 0]'';sigma_w=diag([10e-3 10e-3]);');
disp('>> mu_v=[0 0]'';sigma_v=diag([10e-3 10e-3]);');
disp(' ');
disp('>> x0={mu_x0,sigma_x0};');
disp('>> w={mu_w,sigma_w};');
disp('>> v={mu_v,sigma_v};');
disp(' ');
disp('2) Set the other parameters ');
disp('>> k0 = 0; ');
disp('>> Ts = 0.1; ');
disp('>> TimeUnit=''seconds''; ');
disp(' ');
disp('3) Create the state space model ');
disp('ssObj=ss_DNL_AN(@demo_sim1f,@demo_sim1h,x0,w,v,k0,Ts,TimeUnit);');
disp(' ');
disp('Note that all variables are set explicitly, to make it easier to');
disp('see what is happening. Shorter notations are of course possible.');
mu_x0=[0 5]';
sigma_x0=diag([5 5]);
mu_w=[0 0]';
sigma_w=diag([10e-3 10e-3]);
mu_v=[0 0]';
sigma_v=diag([10e-3 10e-3]);
x0={mu_x0,sigma_x0};
w={mu_w,sigma_w};
v={mu_v,sigma_v};
k0 = 0; 
Ts = 0.1; 
TimeUnit='seconds'; 
ssObj=ss_DNL_AN(@demo_sim1f,@demo_sim1h,x0,w,v,k0,Ts,TimeUnit);
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that the state-space model is ready, we can simulate its states');
disp('and measurements in different ways.');
disp('The fastes way is as follows:');
disp(' ');
disp('>> simObj = sim(ssObj);');
simObj = sim(ssObj);
disp(' ');
disp('This fast calling created the following simulation object:');
disp(' ');
disp('>> simObj');
simObj
disp(' ');
disp('It only contains one sample: one state and one measurement.');
disp('This is the default behaviour if the amount of sample is not');
disp('specified.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now lets try some more samples, for example 10:');
disp(' ');
disp('>> simObj = sim(ssObj,10)');
simObj = sim(ssObj,10)
disp(' ');
disp('The step number index has the following content:');
disp(' ');
disp('>> simObj.k');
simObj.k
disp('It starts from 0 to 9, which are the values of k for each sample.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('If you would like to start from step 5, just change it in the ');
disp('state space model and run the simulation again:');
disp('>> ssObj.k0=5;simObj = sim(ssObj,10);');
disp('>> simObj.k');
ssObj.k0=5;
simObj = sim(ssObj,10);
simObj.k
disp(' ');
disp('Suppose you are not intrested in saving the simulated states. ');
disp('You can change this setting by adding the conf settings: ');
disp('>> simObj = sim(ssObj,10,[],''y'');');
disp('>> simObj');
simObj = sim(ssObj,10,[],'y');
simObj
disp(' ');
disp('Now, the states are not saved. Also note that when calling the method,');
disp('an empty matrix was used to indicate that there are no inputs. ');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Another setting that can be useful is to specify the initial state.');
disp('Until now, the simulation model has simulated the initial state by');
disp('taking a sample from the initial state noise model. For example:');
disp(' ');
disp('>> simObj = sim(ssObj,10);');
disp('>> simObj.x(1,:)');
simObj = sim(ssObj,10);
simObj.x(:,1)
disp(' ');
disp('To specify the inital state manually do the following:');
disp(' ');
disp('>> simObj = sim(ssObj,10,[],'''',[1.4 0]'');');
simObj = sim(ssObj,10,[],'',[1.4 0]');
disp(' ');
disp('Lets see if that was really applied.');
disp(' ');
disp('>> simObj.x(1,:)');
simObj.x(:,1)
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('The last possible setting of the simulation method is the noise');
disp('parameter. When setting this value equal to 0, no noise is added');
disp('during the simulation, when set to 1 only process noise is added,');
disp('when set to 2 only measurement noise is addded and finally, when');
disp('the noise is set to 3 both noise sources are addded. This is of ');
disp('course the default setting.');
disp(' ');
disp('As a final example, lets set the noise to 2 while leaving the other');
disp('at their default values:');
disp(' ');
disp('>> simObj = sim(ssObj,10,[],'''',[]'',2);');
simObj = sim(ssObj,10,[],'',[],2);
disp(' ');
disp('This concludes this demo. ');