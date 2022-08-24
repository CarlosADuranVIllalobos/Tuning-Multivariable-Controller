%% createChannel.m
% Creates a channel structure

% [c] = createChannel(name, measurement_unit, time_unit)
% [c] = createChannel(name, measurement_unit, time_unit, sampling_period, t, y)
%
% name - channel name
% measurement_unit - channel measurement unit
% time_unit - time unit
% t - vector of time instances
% y - vector of measurements
% c has members {name, yUnit, tUnit, t ,y}

function [c] = createChannel(name, measurementUnit, timeUnit, t, y)
switch nargin
    % name and unit
    case 3
        % check if parameters are valid
        if(ischar(name) == 1)
            c.name = name;
        else
            error('Name is not a string!');
        end
        if(ischar(measurementUnit) == 1)
            c.yUnit = measurementUnit;
        else
            error('Measurement unit is not a string!');
        end
        if(ischar(timeUnit) == 1)
            c.tUnit = timeUnit;
        else
            error('Time unit is not a string!');
        end 
    case 5
        % check if parameters are valid
        if(ischar(name) == 1)
            c.name = name;
        else
            error('Name is not a string!');
        end
        if(ischar(measurementUnit) == 1)
            c.yUnit = measurementUnit;
        else
            error('Measurement unit is not a string!');
        end
        if(ischar(timeUnit) == 1)
            c.tUnit = timeUnit;
        else
            error('Time unit is not a string!');
        end 
        if(isreal(t) == 1 && isreal(y) == 1)
            % make t and y column vectors
            [nt, mt] = size(t);
            [ny, my] = size(y);
            if(nt < mt)
                t = t';
                [nt, mt] = size(t);
            end
            if(ny < my)
                y = y';
                [ny, my] = size(y);
            end
            
            % check if t and y have the same number of rows
            if(ny == nt)
                c.t = t(:,1);
                c.y = y(:,1);
            else
                error('t and y dont have the same number of rows!');
            end
        else
            error('t or y are not real arrays!');
        end
        
    otherwise
end
end
%% Copyright
% Stephen Goldrick and Andrei Stefan
% University of Manchester, Newcastle University and Perceptive Engineering
%
% All rights reserved. Copyright (c) Newcastle University, University of Manchester and Perceptive Engineering.

