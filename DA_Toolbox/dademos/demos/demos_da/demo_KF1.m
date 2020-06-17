% Demo about using the KF da technique

clc;

disp('This is a demo of how to apply the Kalman filter KF on a discrete-time ');
disp('state space model to obtain estimated states.');
disp('The example considers a vehicle accelaration and is based on the paper:');
disp(' ''Kalman filtering'' from Dan Simon in Embedded Systems Programming.');
disp(' ');
disp('First a ss_DL state space model is created using the following steps:');
disp(' ');
disp('1) Create the noise models ');
disp('>> Ts=0.1;');
disp('>> mu_x0=[0 0]'';sigma_x0=10^2 * [Ts^4/4 Ts^3/2; Ts^3/2 Ts^2+0.001];');
disp('>> mu_w =[0 0]'';sigma_w= 10^2 * [Ts^4/4 Ts^3/2; Ts^3/2 Ts^2+0.001];');
disp('>> mu_v=0;sigma_v=0.2^2;');
disp('>> x0={mu_x0,sigma_x0};w={mu_w,sigma_w};v={mu_v,sigma_v};');
disp(' ');
disp('2) Set the other parameters ');
disp('>> k0 = 0; samples=600; TimeUnit=''seconds'';');
disp(' ');
disp('3) Set the system matrices and input ');
disp('>> A = [1 Ts; 0 1]; ');
disp('>> B = [Ts^2/2; Ts]; ');
disp('>> C = [1 0]; D = 0;');
disp('>> u = ones(1,samples);');
disp(' ');
disp('4) Create the state space model ');
disp('ssObj=ss_DL(A,B,C,D,x0,w,v,k0,Ts,TimeUnit);');
disp(' ');
disp('Note that all variables are set explicitly, to make it easier to');
disp('follow what is happening. Shorter notations are of course possible.');
Ts=0.1;
mu_x0=[0 0]';sigma_x0=0.2^2 * [Ts^4/4 Ts^3/2; Ts^3/2 Ts^2+0.001];
mu_w =[0 0]';sigma_w= 0.2^2 * [Ts^4/4 Ts^3/2; Ts^3/2 Ts^2+0.001];
mu_v=0;sigma_v=10^2;
x0={mu_x0,sigma_x0};
w={mu_w,sigma_w};
v={mu_v,sigma_v};
k0 = 0; 
samples=600; 
TimeUnit='seconds'; 
A = [1 Ts; 0 1]; 
B = [Ts^2/2; Ts]; 
C = [1 0]; 
D = 0; 
u = ones(1,samples);
ssObj=ss_DL(A,B,C,D,x0,w,v,k0,Ts,TimeUnit);
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that the state-space model is ready, we can simulate its states');
disp('and measurements.');
disp(' ');
disp('>> simObj=sim(ssObj,samples,u);');
simObj=sim(ssObj,samples,u);
disp(' ');
disp('Now that virtual measuremtents are available we can apply the');
disp('da method to estimate the states with KF. Note that the true states');
disp('and the measurements are saved for later usage.');
disp(' ');
disp('>> damObj=da(ssObj,''KF'',simObj,''x y'');');
damObj=da(ssObj,'KF',simObj,'x y');
disp(' ');
disp('Lets see how the results look like. We can easily make a plot of the');
disp('analysis states by using the plot command:');
disp(' ');
disp('>> plot(damObj,''xa'',(1:2));');
plot(damObj,'xa',(1:2));
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('Another possible plot is to compare the true states, analysis states');
disp('and measurements in one plot. To this end, the required numbers of');
disp('measurement values needs to be added since they can differ from the');
disp('state numbers. In this case only the first state is measured so the');
disp('syntax becomes:');
disp(' ');
disp('>> plot(damObj,''xaxy'',1,1);');
plot(damObj,'xaxy',1,1);

disp(' ');
disp('This concludes this demo. ');
disp('press <ENTER> key'); pause
close all;