function xdot = demo_L3f(t,x)
sigma=10;
ro=28;
B=8/3;

        xdot(1,1) = -sigma*(x(1) - x(2));
        xdot(2,1) = x(1)*(ro - x(3)) - x(2);
        xdot(3,1) = x(1)*x(2) - B*x(3);
        
end