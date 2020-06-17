function x1 = demo_VDP_f(x,k,u,~,Ts)
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
alpha=0.2;

      x1(1,1) = x(1) + Ts * x(2);
      x1(2,1) = x(2) + Ts * ( (alpha * (1 - x(1)^2) * x(2) ) - x(1));

end