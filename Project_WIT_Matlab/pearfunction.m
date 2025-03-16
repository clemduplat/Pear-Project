function [x, y] = pearfunction(s)
    
        % Generate x and y coordinates for the custom shape
        x = 0.3 * ((1 + sin(s)) .* cos(s));
        y = -(1 + sin(s)) + 2;
        %x(x < 0) = 0;
        %y = y + 1;
        %x = x/1000;
        %y = y/1000;
end
