% =================================================================== %
% Interactively shows a list of the data assimilation toolbox demos   %
% =================================================================== %
%
% demos_nm_
%
%   demos_nm_gauss
%       - demo_nm_gauss_lti1: Demo about creating and manimulating a linear 
%         time invariant noise model.
%       - demo_nm_gauss_ltv1: Demo about creating and manimulating a linear 
%         time variant noise model.
%       - demo_nm_gauss_handle1: Demo about creating and manimulating a 
%         function handle noise model.
%
% demos_ss_
%
%   demos_ss_D
%       - demo_ss_DL1: Demo about creating and manimulating a linear 
%         time variant state space model.
%       - demo_ss_DNL_AN1: Demo about creating and manimulating a nonlinear 
%         state space model with additive noise.
%       - demo_ss_DNL1: Demo about creating and manimulating a nonlinear 
%         state space model.
%
% demos_sim
%
%   demos_sim
%       - demo_sim1: Demo about using the simulation method to create
%       simulation objects.
%
%   demos_sim_D
%       - demo_sim_D1: Demo about using the methods of the simulation object.
%
% demos_da
%
%       - demo_KF1: Demo about applying the Kalman filter on a vehicle
%       navigation model.
%       - demo_L3: Demo about applying several Ensemble filters on the
%       Lorenz equations.
%       - demo_VDP: Demo about applying the EKF, UKF and EnKF filters on
%       the Van Der Pol oscillator.
%       - demo_TRACK: Demo about applying the EnKF, the UKF and the PF_ASIR
%       filters on a tracking problem with linear state equations and
%       highly nonlinear measurements.
%       - demo_NLSC: Demo about applying the EKF, UKF, EnKF and PF_ASIR on
%       a highly nonlinear scalar model.
%

clc;

disp('To run the demos, enter the name of the demo in the command line.');
disp('(e.g. demo_nm_gauss_lti1)');
disp('The following demos are currently present in the toolbox:');
disp(' ');
disp('	Demos about noise models:');
disp(' ');
disp('      The Gaussian noise models');
disp('      -------------------------');
disp('      - demo_nm_gauss_lti1 -');
disp('      Demo about creating and manimulating a linear time invariant');
disp('      noise model.');
disp(' ');
disp('      - demo_nm_gauss_ltv1 -');
disp('      Demo about creating and manimulating a linear time variant');
disp('      noise model.');
disp(' ');
disp('      - demo_nm_gauss_handle1 -');
disp('      Demo about creating and manimulating a function handle noise');
disp('      model.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('To run the demos, enter the name of the demo in the command line.');
disp('(e.g. demo_ss_DL1)');
disp('The following demos are currently present in the toolbox:');
disp(' ');
disp('	Demos about state space models:');
disp(' ');
disp('      The discrete-time state space models');
disp('      ------------------------------------');
disp('      - demo_ss_DL1 -');
disp('      Demo about creating and manimulating a linear time variant');
disp('      state space model.');
disp(' ');
disp('      - demo_ss_DNL_AN1 -');
disp('      Demo about creating and manimulating a nonlinear state space');
disp('      model with additive noise.');
disp(' ');
disp('      - demo_ss_DNL1 -');
disp('      Demo about creating and manimulating a nonlinear state space');
disp('      model.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('To run the demos, enter the name of the demo in the command line.');
disp('(e.g. demo_ss_DL1)');
disp('The following demos are currently present in the toolbox:');
disp(' ');
disp('	Demos about simulation models:');
disp(' ');
disp('      The discrete-time simulation method');
disp('      ------------------------------------');
disp('      - demo_sim1 -');
disp('      Demo about using the simulation method to create simulation');
disp('      objects.');
disp(' ');
disp('      The discrete-time simulation model');
disp('      ------------------------------------');
disp('      - demo_sim_D1 -');
disp('      Demo about using the methods of the simulation object.');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('To run the demos, enter the name of the demo in the command line.');
disp('(e.g. demo_ss_DL1)');
disp('The following demos are currently present in the toolbox:');
disp(' ');
disp('	Demos about data assimilation methods and models:');
disp(' ');
disp('      The discrete-time data assimilation methods');
disp('      -------------------------------------------');
disp('      - demo_KF1 -');
disp('      Demo about applying the Kalman filter on a vehicle navigation');
disp('      model.');
disp(' ');
disp('      - demo_L3 -');
disp('      Demo about applying several Ensemble filters on the Lorenz ');
disp('      equations.');
disp(' ');
disp('      - demo_VDP -');
disp('      Demo about applying the EKF, UKF and EnKF filters on the Van Der');
disp('      Pol oscillator.');
disp(' ');
disp('      - demo_TRACK -');
disp('      Demo about applying the EnKF, the UKF and the PF_ASIR filters ');
disp('      on a tracking problem with a linear state equation and a highly');
disp('      nonlinear measurement equation.');
disp(' ');
disp('      - demo_NLSC -');
disp('      Demo about applying the EKF, UKF, EnKF and PF_ASIR on a highly');
disp('      nonlinear scalar model.');
disp(' ');

