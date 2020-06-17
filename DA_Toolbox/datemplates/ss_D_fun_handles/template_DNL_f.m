function x1 = template_DNL_f(x,k,u,w,Ts)
%TEMPLATE_DNL_F function handle template for f for model of class ss_DNL
%
%  - Input variable(s) -
%  X: current state vector
%  K: the step number
%  U: input vector
%  W: process noise vector w
%  TS: the sample time
%
%  - Output variable(s) -
%  X1: new state without noise
%
%  - Construction -
%  X1 = TEMPLATE_DNL_F(X,K,U,W,TS) returns the new state X1 with noise w
%  with input U and sample time TS
%

% EXAMPLES:
% linear model
% x1(1,1) = x(1)+2*x(2) + u + w(1);
% x1(2,1) = x(1)+x(2) + w(2);
%
% non-linear model
% x1(1,1) = x(2)*k  + w(1);
% x1(2,1) = -x(1) * x(2)^2 / 2 * x(3)  + w(2);
% x1(3,1) = x(2) + w(3);

end
