classdef unitCell<PLG
    %UNITCELL this class will be a subset of the PLG and will output a unit
    %cell based on a series of methods
    
    properties
        vertices
        connections
        unitType
    end
    
    methods
        function obj = unitCell(unitName,unitType)
            %UNITCELL Constructs this class as loads a unit cell
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    methods (Access=protected)
        %square unit cells
        % unit cells are created by joining several base shapes
        function obj = centreCrossBottom(obj)
            
        end
        
    end
end

