function hJacX = template_DNL_hJacX(x,k,u,v,Ts)
%TEMPLATE_DNL_HJACX function handle template for the Jacobian of h
%                      with respect to X for model of class ss_DNL
%
%  - Input variable(s) -
%  X: current state vector
%  K: the step number
%  U: input vector
%  V: meas. noise vector v
%  TS: the sample time
%
%  - Output variable(s) -
%  HJACX: Jacobian of f with respect to X
%
%  - Construction -
%  HJACX = TEMPLATE_DNL_HJACX(X,K,U,V,TS) returns Jacobian of f with respect
%  to X at step K, with input U and sample time TS
%

% EXAMPLES:
% linear model
% hJacX = [1 0 ; 0 1] ;
%
% non-linear model
% hJacX = [1 k 0];
      

end