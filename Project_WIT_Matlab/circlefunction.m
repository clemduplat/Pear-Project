function [x,y] = circlefunction(bs,s)
% Create a unit circle centered at (0,0) using four segments.
%disp(bs)
%disp(s)
switch nargin
    case 0
        
        x = 4; % four edge segments
        return
    case 1
        
        A = [0,pi/2,pi,3*pi/2; % start parameter values
             pi/2,pi,3*pi/2,2*pi; % end parameter values
             1,1,1,1; % region label to left
             0,0,0,0]; % region label to right
        x = A(:,bs); % return requested columns
        return
    case 2
        x = 0.6*(1+sin(s)).*cos(s);
        y = sin(s) - 1;

        

        x(x < 0) = 0;

        x = x/1000;
        y = y/1000;

end
