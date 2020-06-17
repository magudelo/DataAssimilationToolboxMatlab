function y = template_DNL_h(x,k,u,v,Ts)
%TEMPLATE_DNL_H function handle template for h for model of class ss_DNL
%
%  - Input variable(s) -
%  X: current state vector
%  K: the step number
%  U: input vector
%  V: meas. noise vector v
%  TS: the sample time
%
%  - Output variable(s) -
%  Y: measurement without noise
%
%  - Construction -
%  Y = TEMPLATE_DNL_H(X,K,U,V,TS) returns the new state Y with noise v
%  at step K with input U and sample time TS
%

% EXAMPLES:
% linear model
% y(1,1) = x(1) + v(1);
% y(2,1) = x(2) + v(2);
%
% non-linear model
% y = x(1)+x(2)*k + v;

end
