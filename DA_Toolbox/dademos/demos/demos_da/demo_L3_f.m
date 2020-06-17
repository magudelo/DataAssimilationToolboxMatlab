function x1 = demo_L3_DNL_ANf(x,k,u,~,Ts)
%TEMPLATE_DNL_AN_F function handle template for f for model of class ss_DNL_AN
%
%  - Input variable(s) -
%  X: current state vector
%  K: the step number
%  U: input vector
%  ~: process noise w (does not apply here)
%  TS: the sample time
%
%  - Output variable(s) -
%  X1: new state without noise
%
%  - Construction -
%  X1 = TEMPLATE_DNL_AN_F(X,K,U,~,TS) returns the new state X1 without noise
%  at step K with input U and sample time TS
%

%[TOUT,YOUT] = ode45(ODEFUN,TSPAN,Y0)
[TOUT,x1] = ode45(@demo_L3f,[k*Ts  (k+1)*Ts],x);
x1=x1(end,:)';

end
