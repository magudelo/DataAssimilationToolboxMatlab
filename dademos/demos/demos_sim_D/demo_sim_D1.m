% Demo about using the methods of the simulation object

clc;

disp('This is a demo of how to use the simulation methods of the discrete-');
disp('time simulation object.');
disp('To have a fast start on using the methods we will configure a simulation');
disp('model manually. Lets make an object that contains 200 samples of a ');
disp('state space with 10 states and 5 inputs where only the even states');
disp('are measured as follows');
disp(' ');
disp('1) create 10 states of values 5 to 45 with added noise');
disp(' ');
disp('>> X=randn(10,200); X=X+repmat((0:5:45)'',1,200);');
disp(' ');
X=randn(10,200); X=X+repmat((0:5:45)',1,200);
disp('2) create 5 inputs of values -5 to -20 with added noise');
disp(' ');
disp('>> U=0.1*randn(5,200); U=U-repmat((0:5:20)'',1,200);');
disp(' ');
U=0.1*randn(5,200); U=U-repmat((0:5:20)',1,200);
disp('3) set measurements equal to even state and add noise');
disp(' ');
disp('>> Y=X((2:2:10),:)+0.1*randn(5,200);');
disp(' ');
Y=X((2:2:10),:)+0.1*randn(5,200);
disp('4) create the simulation object with sample time of 2 minutes and an');
disp('   index that starts from 5 till 204');
disp(' ');
disp('>> simObj=sim_D(X,U,Y,(5:204),2,''minutes'')');
disp(' ');
simObj=sim_D(X,U,Y,(5:204),2,'minutes');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('Now that we have a simulation object, we can make plots of it using');
disp('the plot method. The plot method always requires at least three');
disp('input arguments: the simulation object, the configuration string that');
disp('specifies whether the state, input or meas.needs to be shown and');
disp('the specification of which state/input/meas. are requested.');
disp(' ');
disp('As a first example lets plot the states 2-5 and 7-8 using the default ');
disp('plot style (subplots)');
disp(' ');
disp('>> plot(simObj,''x'',[2:5 7:8]) ');
plot(simObj,'x',[2:5 7:8])
disp(' ');
disp('Notice how the plot method automatically arranges the subplots and ');
disp('indicates the different states.');
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('As a second example lets plot all the inputs using the single ');
disp('plot style, which is indicated by adding a 1 in the configuration string');
disp(' ');
disp('>> plot(simObj,''u 1'',(1:5) ');
plot(simObj,'u 1',(1:5))
disp(' ');
disp('Notice how the plot method automatically adds a legend and ');
disp('and a title.');
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('As a third example lets plot all the measurements using the subplot');
disp('plot style while using the time scale for the x-axis by adding a t');
disp('in the configuration string as follows:');
disp(' ');
disp('>> plot(simObj,''y 2 t'',(1:5) ');
plot(simObj,'y 2 t',(1:5))
disp(' ');
disp('Notice how the plot method properly changed the x -axis scale.');
disp(' ');
disp('press <ENTER> key'); pause

clc;close all;
disp('As a final example lets plot the states and measurements using the ');
disp('subplot plot style. Remember that only the even states were measured.');
disp('Hence the command becomes:');
disp(' ');
disp('>> plot(simObj,''xy 2 t'',(2:2:10),(1:5) )');
plot(simObj,'xy 2 t',(2:2:10),(1:5))
disp(' ');
disp('This concludes this demo.');
disp('press <ENTER> key'); pause
close all;
