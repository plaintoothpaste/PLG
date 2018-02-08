classdef unitCell
    %UNITCELL this class will be a subset of the PLG and will output a unit
    %cell based on a series of methods
    
    properties
        vertices
        connections
        unitType
        unitName
        
        diameter
        resolution
        scale
    end
    
    methods
        function obj = unitCell(unitName)
            %UNITCELL Constructs this class as loads a unit cell
            obj.vertices = [];
            obj.connections = [];
            obj.unitType = [];
            obj.unitName = [];
            obj.diameter = [];
            obj.resolution = [];
            
            obj = load(obj,unitName);
        end
        function obj = addUnit(obj,unitName)
            % add a unit cell to the object
            obj = load(obj,unitName);
        end
        function obj = set(obj,varargin)
            % set a value in the PLG that can be edited
            allowable={'diameter','resolution','scale'};
            dataType = {'double','double','double'};
            dataSize = {[1,1],[1,1],[1,3]};
            % setup parser
            p = inputParser();
            for inc = 1:length(allowable)
                name = allowable{inc};
                class = dataType{inc};
                sizer = dataSize{inc};
                f = @(x) validateattributes(x,{class},{'size', sizer,'nonempty'});
                addOptional(p,name,[],f)
            end
            parse(p,varargin{:});
            
            % apply values
            for inc = 1:length(allowable)
                name = allowable{inc};
                result = p.Results.(name);
                if ~isempty(result)
                    obj.(name) = result;
                end
            end
        end
        function obj = beam2facets(obj)
            % converts a beam model to a facet model, primarially done to allow addition of beam and
            % facet unit cells
            % TODO
        end
        function [vertices,connections,transform,name,type] = output(obj)
            % if the model is a beam model then a single strut [0,0,-0.5]->[0,0,0.5] 
            % with its repspective affine transform will be supplied.
            if isempty(obj.diameter)
                error('diameter and resolution must be supplied');
            end
            
            switch obj.unitType
                case 'beam'
                    connections = [1,2];
                    transform = [];
                    vertices = [0,0,-0.5;0,0,0.5];
                    thirdPoint = vertices(1,:)+[1,0,0];
                    forthPoint = vertices(2,:)+[0,1,0];
                    originalPoints = [vertices(1,:),1;vertices(2,:),1;thirdPoint,1;forthPoint,1];
                    name = obj.unitName;
                    type = obj.unitType;
                    obj = scaleBeam(obj);
                    
                    % convert connections to transforms
                    for inc = 1:size(obj.connections,1)
                        point1 = obj.vertices(obj.connections(inc,1),:);
                        point2 = obj.vertices(obj.connections(inc,2),:);
                        vector = point2-point1;
                        u1 = vector/norm(vector);
                        if u1(1)==1 || u1(1)==-1
                            crosser1 = [0,1,0];
                            crosser2 = [0,0,1];
                        elseif u1(2)==1 || u1(2)==-1
                            crosser1 = [0,0,1];
                            crosser2 = [1,0,0];
                        else % perfect z or any other vector
                            crosser1 = [1,0,0];
                            crosser2 = [0,1,0];
                        end
                        v1 = cross(crosser1,u1);
                        v1 = v1/norm(v1);
                        offset = obj.diameter/2*v1;
                        point3 = point1+offset;
                            
                        vector = point1-point2;
                        u1 = vector/norm(vector);
                        v1 = cross(crosser2,u1);
                        v1 = v1/norm(v1);
                        offset = obj.diameter/2*v1;
                        point4 = point2+offset;
                        newPoints = [point1,1;point2,1;point3,1;point4,1;];
                        
                        affine = originalPoints\newPoints;
                        transform = [transform;affine(1,1:3),affine(2,1:3),affine(3,1:3),affine(4,1:3)];
                        
                    end
                case 'facet'
                    transform = [];
            end
                
        end  
    end
    methods (Access=protected)
        function obj = load(obj,unitName)
            % load all xml objects
            fileName = ['unitCell',filesep,unitName,'.xml'];
            xmlObj = xmlread(fileName);
            xmlStructure = unitCell.parseChildNodes(xmlObj);
            
            % convert the xml structure to usefull data
            for inc = 1:length(xmlStructure.Children)
                if strcmp(xmlStructure.Children(inc).Name,'vertices')
                    vertLoc = inc;
                elseif strcmp(xmlStructure.Children(inc).Name,'struts')
                    strutLoc = inc;
                end
            end
            % connections and vertices
            if isempty(obj.vertices)
                strutOffset = 0;
            else
                strutOffset = size(obj.vertices,1);
            end
            verts = xmlStructure.Children(vertLoc).Children;
            for inc = 2:2:length(verts)
                obj.vertices = [obj.vertices;...
                    str2double(verts(inc).Attributes(1).Value),...
                    str2double(verts(inc).Attributes(2).Value),...
                    str2double(verts(inc).Attributes(3).Value)];
            end
            struts = xmlStructure.Children(strutLoc).Children;
            for inc = 2:2:length(struts)
                obj.connections = [obj.connections;...
                    str2double(struts(inc).Attributes(1).Value)+strutOffset,...
                    str2double(struts(inc).Attributes(2).Value)+strutOffset];
            end
            % unit cell type and name
            obj.unitName = [obj.unitName,' ',xmlStructure.Children(strutLoc).Attributes(1).Value];
            if ~isempty(obj.unitType)
                if strcmp(obj.unitType,xmlStructure.Children(strutLoc).Attributes(2).Value);
                    % good load
                else
                    % bad load can not load a beam and facet at the same time
                    error('can not load a beam and a non beam type at the same time add a resolution and diameter then run get facet to allow this')
                end
            else
                % first load
                obj.unitType = xmlStructure.Children(strutLoc).Attributes(2).Value;
            end
        end
        function obj = scaleBeam(obj)
            obj.vertices(:,1) = obj.vertices(:,1)*obj.scale(1);
            obj.vertices(:,2) = obj.vertices(:,2)*obj.scale(2);
            obj.vertices(:,3) = obj.vertices(:,3)*obj.scale(3);
        end
    end
    methods (Static) % xml loading functions
        function children = parseChildNodes(node)
            % Recurse over node children.
            children = [];
            if node.hasChildNodes
                childNodes = node.getChildNodes;
                numChildNodes = childNodes.getLength;
                allocCell = cell(1, numChildNodes);
                
                children = struct(             ...
                    'Name', allocCell, 'Attributes', allocCell,    ...
                    'Data', allocCell, 'Children', allocCell);
                
                for count = 1:numChildNodes
                    theChild = childNodes.item(count-1);
                    children(count) = unitCell.makeStructFromNode(theChild);
                end
            end
        end
        function nodeStruct = makeStructFromNode(node)
            % Create structure of node info.
            
            nodeStruct = struct(                        ...
                'Name', char(node.getNodeName),       ...
                'Attributes', unitCell.parseAttributes(node),  ...
                'Data', '',                              ...
                'Children', unitCell.parseChildNodes(node));
            
            if any(strcmp(methods(node), 'getData'))
                nodeStruct.Data = char(node.getData);
            else
                nodeStruct.Data = '';
            end
        end
        function attributes = parseAttributes(node)
            % Create attributes structure.
            
            attributes = [];
            if node.hasAttributes
                theAttributes = node.getAttributes;
                numAttributes = theAttributes.getLength;
                allocCell = cell(1, numAttributes);
                attributes = struct('Name', allocCell, 'Value', ...
                    allocCell);
                
                for count = 1:numAttributes
                    attrib = theAttributes.item(count-1);
                    attributes(count).Name = char(attrib.getName);
                    attributes(count).Value = char(attrib.getValue);
                end
            end
        end
    end
end

