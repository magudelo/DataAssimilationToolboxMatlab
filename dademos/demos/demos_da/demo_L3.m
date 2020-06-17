% Demo about applying several da techniques on Lorenz equation

clc;

disp('This is a demo of how to apply several data assimilation techniques ');
disp('on the three states Lorenz eqution (L3) to obtain estimated states.');
disp('For more information on the VDP equations see the accompanying thesis.');
disp(' ');
disp('First the ss_DNL_AN state space model is created using the following steps:');
disp(' ');
disp('1) Create the noise models ');
disp('>> x0={diag([1 1 1])};w={diag([16 16 16])};v={diag([4 4 4])};');
disp(' ');
disp('2) Set the other parameters ');
disp('>> k0 = 0; Ts=0.01;TimeUnit=''seconds'';');
disp(' ');
disp('3) Create the state space model ');
disp('ssObj=ss_DNL_AN(@demo_L3_f,@demo_L3_h,x0,w,v,k0,Ts,TimeUnit,...');
disp('          @demo_L3_fjacx,@demo_L3_hjacx);');
disp(' ');
disp('Note that all variables are set explicitly, to make it easier to');
disp('follow what is happening. Shorter notations are of course possible.');
mySeed = 10;
rng(mySeed);             % Set the seed

%x0={diag([1 1 1])};
x0=nm_gauss_lti([-5.8;-5.7;20.5],diag([1 1 1]));

w={diag([1e-6 1e-6 1e-6])};v={diag([2])};
k0 = 0; Ts=0.01;TimeUnit='seconds';
ssObj=ss_DNL_AN(@demo_L3_f,@demo_L3_h,x0,w,v,k0,Ts,TimeUnit,@demo_L3_fjacx,@demo_L3_hjacx);
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that the L3 model is ready, we can simulate its states');
disp('and measurements and save both of them in the sim. object.');
disp(' ');
disp('>> x0init=[3;-3;12];samples=300; conf=''x y'';u=[];');
disp('>> simObj=sim(ssObj,samples,u,conf,x0init);');
%x0init=[3;-3;12];samples=300; conf='x y';u=[];
x0init=[-6;-6;20];samples=600; conf='x y';u=[];

simObj=sim(ssObj,samples,u,conf,x0init);
disp(' ');
disp('We can always quickly check if everything went as expected by making');
disp('a quick plot of the measurements and true states. Lets plot the states');
disp('together with the measurements using the time scale for the x-axis:');
disp(' ');
disp('>> plot(simObj,''xy t'',(1:3),(1:3));');
%plot(simObj,'xy t',(1:3),(1:3));
plot(simObj,'x t',(1:3));
plot(simObj,'y t',(1));

disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('Now that virtual measuremtents are available we can apply the');
disp('da method to estimate the states with DEnKF. Note that x, Pa and y');
disp('are saved for later usage.');
disp(' ');
disp('>> damDEnKF=da(ssObj,''DEnKF'',simObj,''x y Pa'');');
damDEnKF=da(ssObj,'DEnKF',simObj,'x y Pa');
disp(' ');
disp('Lets see how the results look like. We can easily make a plot of the');
disp('true states and the analysis states by using the plot command:');
disp(' ');
disp('>> plot(damDEnKF,''xax'',(1:3));');
plot(damDEnKF,'xax',(1:3));
disp(' ');
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('This looked almost perfect. To have a cleare picture of the errors,');
disp('it might be better to see the errors between the true states and');
disp('the analysis states as follows:');
disp(' ');
disp('>> plot(damDEnKF,''xa-x t'',(1:3));');
plot(damDEnKF,'xa-x t',(1:3));
disp(' ');
disp('Next, to compare the performance of different Ensemble filters');
disp('using the L3 model, the RMSE will be used.');
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('Lets apply the EnKF, the ETKF and the DEnKF ');
disp('with 10 ensemble members, each 3 times on the L3 model to attain');
disp('the average RMSE value as follows:');
disp(' ');
disp('its=3;rmseArray=zeros(3,3); ');
disp('for i=1:its ');
disp('      damEnKF=da(ssObj,''EnKF'',simObj,''x y'',{10}); ');
disp('      damETKF=da(ssObj,''ETKF'',simObj,''x y'',{10}); ');
disp('      damDEnKF=da(ssObj,''DEnKF'',simObj,''x y'',{10}); ');
disp(' ');
disp('      rmseArray(:,1)=rmseArray(:,1)+rmse(damEnKF,(1:3));');
disp('      rmseArray(:,2)=rmseArray(:,2)+rmse(damETKF,(1:3));');
disp('      rmseArray(:,3)=rmseArray(:,3)+rmse(damDEnKF,(1:3));');
disp('end ');
disp('rmseArray=rmseArray./its;rmseArray');
disp(' ');
disp('press <ENTER> key to start iterations. This will take approx. 1 minute');
disp('since the ODE45 solver is slowing down the progress.'); pause

clc;close all;
disp('Please wait...');
its=3;rmseArray=zeros(3,3);
for i=1:its
    
	damEnKF=da(ssObj,'EnKF',simObj,'x y',{10}); 
	damETKF=da(ssObj,'ETKF',simObj,'x y',{10}); 
	damDEnKF=da(ssObj,'DEnKF',simObj,'x y',{10}); 
	rmseArray(:,1)=rmseArray(:,1)+rmse(damEnKF,(1:3));
	rmseArray(:,2)=rmseArray(:,2)+rmse(damETKF,(1:3));
	rmseArray(:,3)=rmseArray(:,3)+rmse(damDEnKF,(1:3));
    disp('iteration:');i
end
rmseArray=rmseArray./its;
rmseArray
disp(' 	EnKF,   	ETKF,    DEnKF');
disp(' ');
disp('As you can see the all algorithms provide very similar estimations');
disp('in this case. ');
disp(' ');
disp('This concludes this demo. ');
disp('press <ENTER> key'); pause
close all;