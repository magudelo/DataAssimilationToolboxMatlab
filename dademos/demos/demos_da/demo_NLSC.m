% Demo about applying several da techniques on a highly nonlinear scalar
% system

clc;

disp('This is a demo of how to apply several data assimilation techniques ');
disp('on a highly nonlinear set of scalar state and measurent eqations');
disp('to obtain estimated states. Although the equations are scalar,');
disp('they are considered as a benchmark for nonlinear filtering due to');
disp('there high nonlinearity.');
disp('For more information about this model see the accompanying thesis.');
disp(' ');
disp('First the ss_DNL_AN state space model is created using the following steps:');
disp(' ');
disp('1) Create the noise models ');
disp('>> x0={0.1,2};w={0,1};v={0,1};');
disp(' ');
disp('2) Set the other parameters ');
disp('>> k0 = 0; Ts=1;TimeUnit=''seconds'';');
disp(' ');
disp('3) Create the state space model ');
disp('ssObj=ss_DNL_AN(@demo_NLSC_f,@demo_NLSC_h,x0,w,v,k0,Ts,TimeUnit,...');
disp('          @demo_NLSC_fjacx,@demo_NLSC_hjacx);');
disp(' ');
disp('Note that all variables are set explicitly, to make it easier to');
disp('follow what is happening. Shorter notations are of course possible.');
mySeed = 10;
rng(mySeed);             % Set the seed
x0={0.1,2};w={0,1};v={0,1};
k0 = 0; Ts=1;TimeUnit='seconds';
ssObj=ss_DNL_AN(@demo_NLSC_f,@demo_NLSC_h,x0,w,v,k0,Ts,TimeUnit,@demo_NLSC_fjacx,@demo_NLSC_hjacx);
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that the L3 model is ready, we can simulate its states');
disp('and measurements and save both of them in the sim. object.');
disp(' ');
disp('>> x0init=0.1;samples=50; conf=''x y'';u=[];');
disp('>> simObj=sim(ssObj,samples,u,conf,x0init);');
x0init=0.1;samples=50; conf='x y';u=[];
simObj=sim(ssObj,samples,u,conf,x0init);
disp(' ');
disp('We can always quickly check if everything went as expected by making');
disp('a quick plot of the measurements and true states. Lets plot the state');
disp('together with the measurement:');
disp(' ');
disp('>> plot(simObj,''xy t'',1,1);');
plot(simObj,'xy t',1,1);
disp(' ');
disp('The highly nonlinear behaviour of both state and measurement are apparent.');
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('Now that virtual measuremtents are available we can apply the');
disp('da method to estimate the states with EKF. Note that x, Pa and y');
disp('are saved for later usage.');
disp(' ');
disp('>> damEKF=da(ssObj,''EKF'',simObj,''x y Pa'');');
damEKF=da(ssObj,'EKF',simObj,'x y Pa');
disp(' ');
disp('Lets see how the results look like. We can easily make a plot of the');
disp('true states and the analysis states with their 95% confidence intervals');
disp('by using the plot command as follows:');
disp(' ');
disp('>> plot(damEKF,''xaxPa'',1);');
plot(damEKF,'xaxPa',1);
disp(' ');
disp('Notice how the estimation fails and in addition the true states are ');
disp('even outside of the confidence intervals at several steps.');
disp('This is an expected behaviour for the EKF when applied on this model.');
disp('Now lets do the same for the EnKF.');
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('For the EnKF, with ensemble size 50, the command is:');
disp(' ');
disp('>> damEnKF=da(ssObj,''EnKF'',simObj,''x y Pa'',{50});');
damEnKF=da(ssObj,'EnKF',simObj,'x y Pa',{50});
disp(' ');
disp('Lets look at the result.');
disp(' ');
disp('>> plot(damEnKF,''xaxPa'',1);');
plot(damEnKF,'xaxPa',1);
disp(' ');
disp('Notice how the estimation has improved a lot and that now the confidence');
disp('intervals are correct.');
disp('Lets see if this is also the case for the PF_ASIR.');
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('For the PF_ASIR, with 50 particles, the command is:');
disp(' ');
disp('>> damASIR=da(ssObj,''PF_ASIR'',simObj,''x y Pa'',{50});');
damASIR=da(ssObj,'PF_ASIR',simObj,'x y Pa',{50});
disp(' ');
disp('Lets look at the result.');
disp(' ');
disp('>> plot(damASIR,''xaxPa'',1);');
plot(damASIR,'xaxPa',1);
disp(' ');
disp('The result looks quite similar as the one of the EnKF, and again the');
disp('confidence intervals are correct.');

disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('Now, lets check the RMSE values of the three algorithms plus the UKF. ');
disp('To obtain a representation of the result, the RMSE experiment is ');
disp('repeated 10 times and only the average result is indicated.');
disp(' ');
disp('its=10;rmseArray=zeros(1,4); ');
disp('for i=1:its ');
disp('      damEKF=da(ssObj,''EKF'',simObj,''x y''); ');
disp('      damEnKF=da(ssObj,''EnKF'',simObj,''x y'',{50}); ');
disp('      damPFASIR=da(ssObj,''PF_ASIR'',simObj,''x y'',{50}); ');
disp('      damUKF=da(ssObj,''UKF'',simObj,''x y''); ');
disp(' ');
disp('      rmseArray(:,1)=rmseArray(:,1)+rmse(damEKF,1);');
disp('      rmseArray(:,2)=rmseArray(:,2)+rmse(damEnKF,1);');
disp('      rmseArray(:,3)=rmseArray(:,3)+rmse(damPFASIR,1);');
disp('      rmseArray(:,4)=rmseArray(:,4)+rmse(damUKF,1);');
disp('end ');
disp('rmseArray=rmseArray./its;rmseArray');
disp(' ');
disp('press <ENTER> key to start iterations.'); pause

clc;close all;
disp('Please wait...');
its=10;rmseArray=zeros(1,4);
for i=1:its
    
	damEKF=da(ssObj,'EKF',simObj,'x y'); 
	damEnKF=da(ssObj,'EnKF',simObj,'x y',{50}); 
	damPFASIR=da(ssObj,'PF_ASIR',simObj,'x y',{50}); 
    damUKF=da(ssObj,'UKF',simObj,'x y'); 
	rmseArray(:,1)=rmseArray(:,1)+rmse(damEKF,1);
	rmseArray(:,2)=rmseArray(:,2)+rmse(damEnKF,1);
	rmseArray(:,3)=rmseArray(:,3)+rmse(damPFASIR,1);
    rmseArray(:,4)=rmseArray(:,4)+rmse(damUKF,1);
    disp('iteration:');i
end
rmseArray=rmseArray./its;
rmseArray
disp(' 	EKF     	EnKF     PF_ASIR    UKF');
disp(' ');
disp('As expected, the EKF has the bad performance. The UKF has the worst');
disp('performance in this cas. Still, remind that it only uses 3 sigma points.');
disp('The particle filter provides the best performance in this case. ');
disp(' ');
disp('This concludes this demo. ');
close all;