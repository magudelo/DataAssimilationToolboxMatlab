function fJacW = template_DNL_fJacW(x,k,u,w,Ts)
%TEMPLATE_DNL_FJACW function handle template for the Jacobian of f
%                      with respect to W for model of class ss_DNL
%
%  - Input variable(s) -
%  X: current state vector
%  K: the step number
%  U: input vector
%  W: process noise vector w
%  TS: the sample time
%
%  - Output variable(s) -
%  FJACW: Jacobian of f with respect to W
%
%  - Construction -
%  FJACW = TEMPLATE_DNL_FJACW(X,K,U,W,TS) returns Jacobian of f with respect
%  to W at step K, with input U and sample time TS
%

% EXAMPLES:
% linear model
% fJacW = [1 0 ; 0 1] ;
%
% non-linear model
% fJacW = [1 0 0;0 1 0;0 0 1];
      

end