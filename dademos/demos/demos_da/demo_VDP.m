% Demo about applying several da techniques on Van der Pol Oscillator

clc;

disp('This is a demo of how to apply several data assimilation techniques ');
disp('on the Van Der Pol Oscillator (VDP) to obtain estimated states.');
disp('For more information on the VDP equations see the accompanying thesis.');
disp(' ');
disp('First the ss_DNL_AN state space model is created using the following steps:');
disp(' ');
disp('1) Create the noise models ');
disp('>> mu_x0=[0 5]'';sigma_x0=diag([5, 5]);');
disp('>> mu_w=[0 0]'';sigma_w=diag([10e-3 10e-3]);');
disp('>> mu_v=[0 0]'';sigma_v=diag([10e-3 10e-3]);');
disp('>> x0={mu_x0,sigma_x0};w={mu_w,sigma_w};v={mu_v,sigma_v};');
disp(' ');
disp('2) Set the other parameters ');
disp('>> k0 = 0; Ts=0.1;TimeUnit=''seconds'';');
disp(' ');
disp('3) Create the state space model ');
disp('ssObj=ss_DNL_AN(@demo_VDP_f,@demo_VDP_h,x0,w,v,k0,Ts,TimeUnit,...');
disp('          @demo_VDP_fjacx,@demo_VDP_hjacx);');
disp(' ');
disp('Note that all variables are set explicitly, to make it easier to');
disp('follow what is happening. Shorter notations are of course possible.');
mySeed = 10;
rng(mySeed);             % Set the seed
mu_x0=[0 5]';sigma_x0=diag([5, 5]);
mu_w=[0 0]';sigma_w=diag([10e-2 10e-2]);
mu_v=[0 0]';sigma_v=diag([10e-2 10e-2]);
x0={mu_x0,sigma_x0};w={mu_w,sigma_w};v={mu_v,sigma_v};
k0 = 0; Ts=0.1;TimeUnit='seconds';
ssObj=ss_DNL_AN(@demo_VDP_f,@demo_VDP_h,x0,w,v,k0,Ts,TimeUnit,@demo_VDP_fjacx,@demo_VDP_hjacx);
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that the VDP model is ready, we can simulate its states');
disp('and measurements and save both of them in the sim. object.');
disp(' ');
disp('>> x0init=[1.4 0]'';samples=500; conf=''x y'';u=[];');
disp('>> simObj=sim(ssObj,samples,u,conf,x0init);');
x0init=[1.4 0]';samples=500; conf='x y';u=[];
simObj=sim(ssObj,samples,u,conf,x0init);
disp(' ');
disp('We can always quickly check if everything went as expected by making');
disp('a quick plot of the measurements and true states as follows:');
disp(' ');
disp('>> plot(simObj,''xy'',(1:2),(1:2));');
plot(simObj,'xy',(1:2),(1:2));
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
disp('true states and the analysis states by using the plot command:');
disp(' ');
disp('>> plot(damEKF,''xax'',(1:2));');
plot(damEKF,'xax',(1:2));
disp(' ');
disp('This looks almost perfect. To compare the performance of diffent filters');
disp('in this case, the RMSE will be used.');
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('Lets apply the EKF, the UKF (has 5 sigma points) and the EnKF ');
disp('with 10 ensemble members each 5 times on the VDP model to attain');
disp('the average RMSE value as follows:');
disp(' ');
disp('its=5;rmseArray=zeros(2,3); ');
disp('for i=1:its ');
disp('      damEKF=da(ssObj,''EKF'',simObj,''x y Pa''); ');
disp('      damUKF=da(ssObj,''UKF'',simObj,''x y Pa''); ');
disp('      damEnKF=da(ssObj,''EnKF'',simObj,''x y Pa'',{10}); ');
disp(' ');
disp('      rmseArray(:,1)=rmseArray(:,1)+rmse(damEKF,(1:2));');
disp('      rmseArray(:,2)=rmseArray(:,2)+rmse(damUKF,(1:2));');
disp('      rmseArray(:,3)=rmseArray(:,3)+rmse(damEnKF,(1:2));');
disp('end ');
disp('rmseArray=rmseArray./its;rmseArray');
disp(' ');
disp('press <ENTER> key to start iterations'); pause

clc;close all;
its=5;rmseArray=zeros(2,3);
for i=1:its
    
    damEKF=da(ssObj,'EKF',simObj,'x y Pa');
    damUKF=da(ssObj,'UKF',simObj,'x y Pa');
    damEnKF=da(ssObj,'EnKF',simObj,'x y Pa',{10});
    
    rmseArray(:,1)=rmseArray(:,1)+rmse(damEKF,(1:2));
    rmseArray(:,2)=rmseArray(:,2)+rmse(damUKF,(1:2));
    rmseArray(:,3)=rmseArray(:,3)+rmse(damEnKF,(1:2));
    disp('iteration:');i
end
rmseArray=rmseArray./its;
rmseArray
disp(' 	EKF,        UKF,    EnKF');
disp(' ');
disp('As you can see the UKF algorithm provides the best estimation in this ');
disp('case. Of course, only when increasing the ensemble size of the EnKF');
disp('its estimation would improve.');
disp(' ');
disp('This concludes this demo. ');
disp('press <ENTER> key'); pause
close all;