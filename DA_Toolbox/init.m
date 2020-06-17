function init(save)
%INIT adds the data assimilation toolbox paths to the Matlab search path
%
%  - Input variable(s) -
%  SAVE: integer value to indicate whether the paths need to be saved for
%  future sessions. (0=do not save)
%
%  - Construction -          
%  INIT(SAVE) adds the data assimilation toolbox paths to the Matlab search 
%  path and does not save it when SAVE=0.
%
%  INIT() adds and saves the data assimilation toolbox paths to the Matlab 
%  search path.

    if nargin == 0
        save = 1;
    end
        
        % add the different folders to the search path
        currDir=pwd;
        addpath(currDir);
        addpath(genpath(strcat(currDir,'\da')));
        addpath(genpath(strcat(currDir,'\dademos')));
        addpath(genpath(strcat(currDir,'\daguis')));
        addpath(genpath(strcat(currDir,'\dahelp')));
        addpath(genpath(strcat(currDir,'\daresource')));
        addpath(genpath(strcat(currDir,'\dautils')));
        addpath(genpath(strcat(currDir,'\datemplates')));        

    if save ~= 0
        savepath
    end
    
    clc;

    disp('The Data Assimilation Toolbox has been properly installed. ');
    disp('  ');
    disp('When first using the toolbox it is recommended to run through');
    disp('the included demos. To see the complete list of the available');
    disp('demos enter ''dademos'' in the command line.');
    disp(' ');
    disp('Also the manual of the data assimilation toolbox which can be');
    disp('found under the folder ''dahelp'' is very instructive.');
    disp('Additional information can be found in the help files of the');    
    disp('classes and functions. A list of the most important help files');    
    disp('can be retrieved by entering ''dahelp'' in the command line.');    
    disp(' ');
    disp('Finally, several templates of functions are available under the');    
    disp('folder datemplates.');    
    disp(' ');
    disp('We are confident that you will enjoy using the DA Toolbox.');
   
end