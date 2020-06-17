function sigma = template_sigma(k,Ts)
%TEMPLATE_SIGMA  function handle template to retrieve sigma for class nm_gauss_handle
%
%  - Input variable(s) -
%  K: the step number
%  TS: the sample time
%
%  - Output variable(s) -
%  SIGMA: the calculated covariance matrix  
%
%  - Construction -
%  SIGMA = TEMPLATE_SIGMA(K,TS) returns the covariance matrix SIGMA of the 
%  noise model with sample time TS at step k
%

% This function can be configured as desired, but keep its syntax!!
%
% EXAMPLES: 
% 
% a bivariate covariance matrix whos variances depend on the sampling time
% sigma =[Ts 0 ; 0 Ts];
%
% shorter syntax of previous example
% sigma =[Ts ; Ts];
%
% mixture with step number
% sigma =[Ts+0.01*k ; Ts];

end