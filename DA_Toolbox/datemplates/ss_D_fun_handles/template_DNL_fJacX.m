function fJacX = template_DNL_fJacX(x,k,u,w,Ts)
%TEMPLATE_DNL_FJACX function handle template for the Jacobian of f
%                      with respect to X for model of class ss_DNL
%
%  - Input variable(s) -
%  X: current state vector
%  K: the step number
%  U: input vector
%  W: process noise vector w
%  TS: the sample time
%
%  - Output variable(s) -
%  FJACX: Jacobian of f with respect to X
%
%  - Construction -
%  FJACX = TEMPLATE_DNL_FJACX(X,K,U,W,TS) returns Jacobian of f with respect
%  to X at step K, with input U and sample time TS
%

% EXAMPLES:
% linear model
% fJacX = [1 2 ; 1 1] ;
%
% non-linear model
% fJacX = [0 k 0;-1 x(2) 1;0 0 0];
      

end