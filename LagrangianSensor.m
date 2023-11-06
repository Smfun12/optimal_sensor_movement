classdef LagrangianSensor
    properties
        x
        y
        path
        color
        idx
        idxy
    end
    
    methods
        % Constructor method
        function obj = LagrangianSensor(x, y, path, color, idx, idxy)
            if nargin > 0
                obj.x = x;
                obj.y = y;
                obj.path = path;
                obj.color = color;
                obj.idx = idx;
                obj.idxy = idxy;
            end
        end
    end
end