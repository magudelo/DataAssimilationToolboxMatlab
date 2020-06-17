function x1 = demo_DNL_AN_f(x,k,u,~,Ts)

      x1(1,1) = x(1)+2*x(2) + u;
      x1(2,1) = x(1)+x(2);

end
