function [string] = ns(number,precision)
% Shorthand way for doing num2str

if nargin<2
    [string] = num2str(number);
else
    [string] = num2str(number,precision);
end

end
