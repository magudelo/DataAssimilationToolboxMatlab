function [rows,cols]=numSubplots(n)
%NUMSUBPLTS determines a proper size for subplots in a nonfancy way
%
%  - Input variable(s) -
%  N: required amount of subplots
%
%  - Output variable(s) -
%  ROWS:amount of rows
%  COLS:amount of columns
%
%  - Construction -
%  [ROWS, COLS] = NUMSUBPLTS(N) determines a proper size for subplots in a 
%  nonfancy way. Current maximum =20. (Should be more than enough)

    switch n
        case 1; rows=1;cols=1;
        case 2; rows=1;cols=2;
        case 3; rows=3;cols=1;
        case 4; rows=2;cols=2;
        case {5,6}; rows=2;cols=3;
        case 7,8; rows=2;cols=4;
        case 9; rows=3;cols=3;
        case {10,11,12}; rows=3;cols=4;
        case {13,14,15}; rows=3;cols=5;
        case {16,17,18,19,20}; rows=4;cols=5;
        otherwise; rows=0;cols=0;
    end
    
end


