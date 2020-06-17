function y = demo_L3_DNL_ANh(x,k,u,~,Ts)
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

C = [1 0 0
     0 1 0
     0 0 1];
C=[1 0 0];
    
 
y=C*x;
 
end
