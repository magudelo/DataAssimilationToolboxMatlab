% Demo about applying several da techniques on a tracking problem

clc;

disp('This is a demo of how to apply several data assimilation techniques ');
disp('on the four states tracking problem to obtain estimated states.');
disp('The tracking problem is interesting for its highly nonlinear measurements.');
disp('In fact two variants are taken under consideration: one measurement');
disp('with a radar which yields distance and angle and one using triangulation');
disp('which provides two different distances. The four states are respectively');
disp('the x position, the y position and the velocities vx and vy. The state');
disp('equation for this model is identical for both radar and triangulation.');
disp('It is important to note that the state equation is linear.');
disp('For more information on the equations see the accompanying thesis.');
disp(' ');
disp('First the ss_DNL_AN state space model is created using the following steps:');
disp(' ');
disp('1) Create the noise models ');
disp('>> x0={[-200;200;4;0]};');
disp('>> w={diag([0.0000001;0.0000001;0.5;0.5])};');
disp('>> vR={diag([200;0.003])};vT={diag([200;200])};');
disp(' ');
disp('2) Set the other parameters ');
disp('>> k0 = 0; Ts=1;TimeUnit=''seconds'';');
disp(' ');
disp('3) Create the state space model with radar measurement');
disp('ssObjR=ss_DNL_AN(@demo_TRACK_f,@demo_TRACK_R_h,x0,w,vR,k0,Ts,TimeUnit');
disp('4) Create the state space model with triangulation measurement');
disp('ssObjT=ss_DNL_AN(@demo_TRACK_f,@demo_TRACK_T_h,x0,w,vT,k0,Ts,TimeUnit');
disp(' ');
disp('Note that all variables are set explicitly, to make it easier to');
disp('follow what is happening. Shorter notations are of course possible.');
mySeed = 10;
rng(mySeed);             % Set the seed
x0={[-200;200;4;0]};w={diag([0.0000001;0.0000001;0.5;0.5])};
vR={diag([200;0.003])};vT={diag([200;200])};
k0 = 0; Ts=0.01;TimeUnit='seconds';
ssObjR=ss_DNL_AN(@demo_TRACK_f,@demo_TRACK_R_h,x0,w,vR,k0,Ts,TimeUnit);
ssObjT=ss_DNL_AN(@demo_TRACK_f,@demo_TRACK_T_h,x0,w,vT,k0,Ts,TimeUnit);
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that the tracking models are ready, we can simulate their states');
disp('and measurements and save both of them in the sim. object.');
disp(' ');
disp('>> x0init=[-200;200;4;0];samples=50; conf=''x y'';u=[];');
disp('>> simObjR=sim(ssObjR,samples,u,conf,x0init);');
disp('>> simObjT=sim(ssObjT,samples,u,conf,x0init);');
x0init=[-200;200;4;0];samples=50; conf='x y';u=[];
simObjR=sim(ssObjR,samples,u,conf,x0init);
simObjT=sim(ssObjT,samples,u,conf,x0init);
disp(' ');
disp('We can  check if everything went as expected by making a quick');
disp('plot of the true states of either one of the models:');
disp(' ');
disp('>> plot(simObjR,''x t'',(1:4));');
plot(simObjR,'x t',(1:4));
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('Now that virtual measuremtents are available we can apply the');
disp('da method to estimate the states with the UKF. Note that x and y');
disp('are saved for later usage. So, for the radar model:');
disp(' ');
disp('>> damRUKF=da(ssObjR,''UKF'',simObjR,''x y'');');
damRUKF=da(ssObjR,'UKF',simObjR,'x y');
disp(' ');
disp('Lets see how the results look like. We can easily make a plot of the');
disp('true states and the analysis states by using the plot command:');
disp(' ');
disp('>> plot(damRUKF,''xax'',(1:4));');
plot(damRUKF,'xax',(1:4));
disp(' ');
disp('It looks like it is difficult to get good estimations. Lets check how');
disp('the estimation is for the triangulation model.');
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('Now that virtual measuremtents are available we can apply the');
disp('da method to estimate the states with the UKF. Note that x and y');
disp('are saved for later usage. So, for the radar model:');
disp(' ');
disp('>> damTUKF=da(ssObjT,''UKF'',simObjT,''x y'');');
damTUKF=da(ssObjT,'UKF',simObjT,'x y');
disp(' ');
disp('Lets see how the results look like. We can easily make a plot of the');
disp('true states and the analysis states by using the plot command:');
disp(' ');
disp('>> plot(damTUKF,''xax'',(1:4));');
plot(damTUKF,'xax',(1:4));
disp(' ');
disp('Still not perfect, but it seems better. Lets now check the accuracy of');
disp('other possible filters using the RMSE.');
disp(' ');
disp('press <ENTER> key'); pause
%%
clc;close all;
disp('The following filters are compared for both the radar as the triangulation');
disp('model: UKF, EnKF with 50 ensemble members and PF_ASIR with 200 particles.');
disp(' ');
disp('Lets start with the radar problem. The average RMSE value of three ');
disp('experiments is checked:');
disp(' ');
disp('its=3;rmseArray=zeros(3,3); ');
disp('for i=1:its ');
disp('      damUKF=da(ssObjR,''UKF'',simObjR,''x y''); ');
disp('      damEnKF=da(ssObjR,''EnKF'',simObjR,''x y'',{50}); ');
disp('      damPFASIR=da(ssObjR,''PF_ASIR'',simObjR,''x y'',{200}); ');
disp(' ');
disp('      rmseArray(:,1)=rmseArray(:,1)+rmse(damUKF,(1:3));');
disp('      rmseArray(:,2)=rmseArray(:,2)+rmse(damEnKF,(1:3));');
disp('      rmseArray(:,3)=rmseArray(:,3)+rmse(damPFASIR,(1:3));');
disp('end ');
disp('rmseArray=rmseArray./its;rmseArray');
disp(' ');
disp('press <ENTER> key to start iterations.');pause

clc;close all;
disp('Please wait...');
its=3;rmseArray=zeros(3,3);
for i=1:its
    
	damUKF=da(ssObjR,'UKF',simObjR,'x y'); 
	damEnKF=da(ssObjR,'EnKF',simObjR,'x y',{50}); 
	damPFASIR=da(ssObjR,'PF_ASIR',simObjR,'x y',{200}); 
	rmseArray(:,1)=rmseArray(:,1)+rmse(damUKF,(1:3));
	rmseArray(:,2)=rmseArray(:,2)+rmse(damEnKF,(1:3));
	rmseArray(:,3)=rmseArray(:,3)+rmse(damPFASIR,(1:3));
    disp('iteration:');i
end
rmseArray=rmseArray./its;
rmseArray
disp(' 	UKF,   	EnKF,    PF_ASIR');
disp(' ');
disp('As you can see the UKF algorithm provides the best estimations with,');
disp('only 9 sigma points which are selected deterministicly. The EnKF');
disp('performs slightly less good but needs 50 ensemble members. The ');
disp('particle filter clearly performs the worst, even with 200 Particles.');
disp('Note however, that the UKF is not suited for large-scale problems.');
disp(' ');
disp('Lets try the same setup for the traingulation model. ');
disp(' ');
disp('press <ENTER> key to start iterations.');pause

clc;close all;
disp('Please wait...');
its=3;rmseArray=zeros(3,3);
for i=1:its
    
	damUKF=da(ssObjT,'UKF',simObjT,'x y'); 
	damEnKF=da(ssObjT,'EnKF',simObjT,'x y',{50}); 
	damPFASIR=da(ssObjT,'PF_ASIR',simObjT,'x y',{200}); 
	rmseArray(:,1)=rmseArray(:,1)+rmse(damUKF,(1:3));
	rmseArray(:,2)=rmseArray(:,2)+rmse(damEnKF,(1:3));
	rmseArray(:,3)=rmseArray(:,3)+rmse(damPFASIR,(1:3));
    disp('iteration:');i
end
rmseArray=rmseArray./its;
rmseArray
disp(' 	UKF,   	EnKF,    PF_ASIR');
disp(' ');
disp('For this model, the estimations of the algorithms are similar.');
disp(' ');
disp('This concludes this demo. ');
disp('press <ENTER> key'); pause
close all;