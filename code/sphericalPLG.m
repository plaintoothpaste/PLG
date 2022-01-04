classdef sphericalPLG < PLG
    % sphericalPLG add cartesion to 3D polar transform
    
    properties
        % see PLG
    end
    
    methods
        function obj = cart2polar(obj)
            % cart2polar Converts the spacial system
            %    Returns a modified PLG object.
            rad = obj.vertices(:,1);
            theta = obj.vertices(:,2);
            omega = obj.vertices(:,3);
            
            obj.vertices(:,1) = rad.*cos(theta).*sin(omega);
            obj.vertices(:,2) = rad.*sin(theta).*sin(omega);
            obj.vertices(:,3) = rad.*cos(omega);
        end
    end
end