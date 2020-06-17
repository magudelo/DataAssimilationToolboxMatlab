function value = iszeromean(x)
%ISZEROMEAN Check if vector is zero mean (all zeros)
%
%  - Input variable(s) -
%  X: vector
%
%  - Output variable(s) - 
%  VALUE: Indicates if X is zero mean (1=zero mean, 0=not zero mean)
%
%  - Construction -
%  VALUE = ISZEROMEAN(X) indicates if X is zero mean.
%

value = ~any(x);