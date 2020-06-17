function str=f2str(func)
%F2STR Function handle to string conversion 
%
%  - Input variable(s) -
%  FUNC: a function handle
%
%  - Output variable(s) -
%  STR: sting containing the function handle
%
%  - Construction -
%  STR = F2STR(FUNC) returns a sting containing the function handle FUNC. 

    if isempty(func)
        str='- empty function -';
    elseif isequal(func,0)
        str='- zero function -';
    else
        str = func2str(func);
        if str(1)~='@'
            str=['@',str];
        end
    end

end