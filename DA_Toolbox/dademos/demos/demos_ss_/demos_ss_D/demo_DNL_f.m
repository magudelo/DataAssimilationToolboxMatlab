function x1 = demo_DNL_f(x,k,u,w,Ts)

      x1(1,1) = x(1)+2*x(2) + u + w(1);
      x1(2,1) = x(1)+x(2)  + w(2);

end
