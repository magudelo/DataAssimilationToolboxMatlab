function hJacV = template_DNL_hJacV(x,k,u,v,Ts)
%TEMPLATE_DNL_HJACX function handle template for the Jacobian of h
%                      with respect to V for model of class ss_DNL
%
%  - Input variable(s) -
%  X: current state vector
%  K: the step number
%  U: input vector
%  V: meas. noise vector v
%  TS: the sample time
%
%  - Output variable(s) -
%  HJACV: Jacobian of f with respect to V
%
%  - Construction -
%  HJACV = TEMPLATE_DNL_HJACV(X,K,U,V,TS) returns Jacobian of f with respect
%  to V at step K, with input U and sample time TS
%

% EXAMPLES:
% linear model
% hJacV = [1 0 ; 0 1] ;
%
% non-linear model
% hJacV = 1;
      

end