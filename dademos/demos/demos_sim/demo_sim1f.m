function x1 = vdp1_f(x,t,u,~,Ts)

alpha=0.2;

      x1(1,1) = x(1) + Ts * x(2);
      x1(2,1) = x(2) + Ts * ( (alpha * (1 - x(1)^2) * x(2) ) - x(1));

end