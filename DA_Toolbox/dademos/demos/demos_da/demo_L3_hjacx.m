function hJacX = demo_L3_DNL_ANhJacX(x,k,u,~,Ts)
%TEMPLATE_DNL_AN_HJACX function handle template for the Jacobian of h
%                      with respect to X for model of class ss_DNL_AN
%
%  - Input variable(s) -
%  X: current state vector
%  K: the step number
%  U: input vector
%  ~: meas. noise v (does not apply here)
%  TS: the sample time
%
%  - Output variable(s) -
%  HJACX: Jacobian of f with respect to X
%
%  - Construction -
%  HJACX = TEMPLATE_DNL_AN_HJACX(X,K,U,~,TS) returns Jacobian of f with respect
%  to X at step K, with input U and sample time TS
%

hJacX = [1 0 0
         0 1 0
         0 0 1];
      

end