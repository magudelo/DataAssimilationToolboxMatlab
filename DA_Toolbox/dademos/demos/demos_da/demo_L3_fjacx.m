function fJacX = demo_fbodyfJacX(x,k,u,~,Ts)
%TEMPLATE_DNL_AN_FJACX function handle template for the Jacobian of f
%                      with respect to X for model of class ss_DNL_AN
%
%  - Input variable(s) -
%  X: current state vector
%  K: the step number
%  U: input vector
%  ~: process noise w (does not apply here)
%  TS: the sample time
%
%  - Output variable(s) -
%  FJACX: Jacobian of f with respect to X
%
%  - Construction -
%  FJACX = TEMPLATE_DNL_AN_FJACX(X,K,U,~,TS) returns Jacobian of f with respect
%  to X at step K, with input U and sample time TS
%
sigma=10;
ro=28;
B=8/3;

fJacX = [-sigma     0       0
         ro-x(3)    -1      x(1)
         x(2)       x(1)    -B];

end