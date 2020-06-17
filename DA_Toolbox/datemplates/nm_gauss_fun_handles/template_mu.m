function mu = template_mu(k,Ts)
%TEMPLATE_MU  function handle template to retrieve mu for class nm_gauss_handle
%
%  - Input variable(s) -
%  K: the step number
%  TS: the sample time
%
%  - Output variable(s) -
%  MU: the calculated mean  
%
%  - Construction -
%  MU = TEMPLATE_MU(K,TS) returns the mean MU of the noise model with
%  sample time TS at step k
%

% This function can be configured as desired, but keep its syntax!!
%
% EXAMPLES: 
% 
% a bivariate mean where the second variable increases proportionaly with k
% mu=[0 ;0.01*k*Ts];
%
% a LTI bivariate zero mean
% mu=[0 ; 0];

end