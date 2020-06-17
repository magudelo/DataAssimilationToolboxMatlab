% ============================================================= %
% Interactively shows a list of the most important help files   %
% ============================================================= %
%
% Noise models:
%       - nm_gauss_lti1: properties, construction and methods
%       - nm_gauss_ltv1: properties, construction and methods
%       - nm_gauss_handle1: properties, construction and methods
%
% State space models:
%
%       - ss_D: da and sim methods
%       - ss_DL: properties, construction and methods
%       - ss_DNL_AN: properties, construction and methods
%       - ss_DNL: properties, construction and methods
%
% Simulation models:
%
%       - sim_D: properties, construction and methods
%
% Data assimilation models:
%
%       - dam_D: properties, construction and methods
%

clc;

disp('To see the help files, enter ''help'' and the name of the file in the ');
disp('command line (e.g. help nm_gauss_lti).');
disp('All the help files explain the different properties, construction techniques');
disp('and methods of the corresponding class. To have more detailed information');
disp('about a method, click on its highlighted name. The help files of methods,');
disp('or functions in general, explain the in- and output parameters of the ');
disp('function and the syntax to call it.');
disp(' ');
disp('	Help files for noise models:');
disp(' ');
disp('      The Gaussian noise models');
disp('      -------------------------');
disp('      - nm_gauss_lti - The linear time invariant Gaussian noise model');
disp('      - nm_gauss_ltv - The linear time variant Gaussian noise model');
disp('      - nm_gauss_handle - The Gaussian noise model with function handles');
disp(' ');
disp('	Help files for state space models:');
disp(' ');
disp('      The discrete-time state space models');
disp('      ------------------------------------');
disp('      - ss_DL - The linear discrete-time state space model');
disp('      - ss_DNL_AN - The nonlinear discrete-time state space model');
disp('                    with additive noise');
disp('      - ss_DNL - The nonlinear discrete-time state space model');
disp(' ');
disp('      - ss_D - The superclass of all discrete-time state space models.');
disp('               Its methods are the simulation method and all data ');
disp('               assimilation techniques. ');
disp(' ');
disp('press <ENTER> key'); pause

clc;
disp('To see the help files, enter ''help'' and the name of the file in the ');
disp('command line (e.g. help nm_gauss_lti).');
disp('All the help files explain the different properties, construction techniques');
disp('and methods of the corresponding class. To have more detailed information');
disp('about a method, click on its highlighted name. The help files of methods,');
disp('or functions in general, explain the in- and output parameters of the ');
disp('function and the syntax to call it.');
disp(' ');
disp('	Help files for simultation models:');
disp(' ');
disp('      The discrete-time simulation models');
disp('      -----------------------------------');
disp('      - sim_D - The discrete-time simulation model.');
disp(' ');
disp('	Help files for data assimilation models:');
disp(' ');
disp('      The discrete-time data assimilation models');
disp('      ------------------------------------------');
disp('      - dam_D - The discrete-timecdata assimilation model.');
disp(' ');

