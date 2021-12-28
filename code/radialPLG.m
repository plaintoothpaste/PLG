classdef radialPLG < PLG
    %radialPLG add cartesion to radial transform
    
    properties
        % see PLG
    end
    
    methods
        function obj = cart2radial(obj)
            % cart2radial Converts the spacial system
            %    Returns a modified PLG object.
            rad = obj.vertices(:,1);
            theta = obj.vertices(:,2);
          
            obj.vertices(:,1) = rad.*cos(theta);
            obj.vertices(:,2) = rad.*sin(theta);
        end
    end
end