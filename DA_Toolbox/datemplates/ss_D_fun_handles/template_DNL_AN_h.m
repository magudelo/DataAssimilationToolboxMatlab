function y = template_DNL_AN_h(x,k,u,~,Ts)
%TEMPLATE_DNL_AN_H function handle template for h for model of class ss_DNL_AN
%
%  - Input variable(s) -
%  X: current state vector
%  K: the step number
%  U: input vector
%  ~: meas. noise v (does not apply here)
%  TS: the sample time
%
%  - Output variable(s) -
%  Y: measurement without noise
%
%  - Construction -
%  Y = TEMPLATE_DNL_AN_H(X,K,U,~,TS) returns the new state Y without noise 
%  at step K with input U and sample time TS
%

% EXAMPLES:
% linear model
% y(1,1) = x(1);
% y(2,1) = x(2);
%
% non-linear model
% y = x(1)+x(2)*k;

end
