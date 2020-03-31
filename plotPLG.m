classdef plotPLG
    % The plotPLG program will plot custom files with desired settings
    properties (SetAccess=protected)
        f
        a
        c
        cMap % colourMap
        gaStruts % a grapics array for the struts
        scVerts % a scatter3 of the verts
    end
    properties (SetAccess=private)
        vertices
        struts
        strutDiameter
        sphereDiameter
    end
    methods
        function obj = plotPLG(fName)
            % load the data and prep a basic plot
            obj = load(obj,fName);
            obj = figureSetupAndPosition(obj);
            obj = defaultPlot(obj);
        end
    end
    methods (Access=private) % called when plotPLG is used
        function obj = load(obj,file)
            % load a custom beam input file to generate a lattice structure
            [~,~,ext] = fileparts(file);
            switch ext
                case '.custom'
                    % standard input
                    data = csvread(file);
                otherwise
                    error('not a suitable load format');
            end
            numNodes=data(1,1);
            numLinks=data(2,1);
            % data
            obj.vertices = data(3:numNodes+2,1:3);
            obj.struts    = data(numNodes+3:numNodes+numLinks+2,1:2);
            obj.strutDiameter = data(numNodes+3:numNodes+numLinks+2,3);
            if size(data,2)==3
                % no sphere diameter supplied
                obj.sphereDiameter = zeros(numNodes,1);
            else
                % sphere diameter supplied
                obj.sphereDiameter = data(3:numNodes+2,4);
            end
        end
        function obj = figureSetupAndPosition(obj)
            % create the figure (f) and axes (a)
            ranger = [min(obj.vertices);max(obj.vertices)];
            d = [obj.strutDiameter;obj.sphereDiameter];
            rangerD = [min(d), max(d)];
            if rangerD(1)==rangerD(2)% in case value 1 and 2 match
                rangerD(1) = rangerD(1)-0.1; 
                rangerD(2) = rangerD(2)-0.1;
            end
            
            obj.f = figure;
            obj.f.Units	= 'centimeters';
            obj.f.Position = [3,3,12,12];
            obj.f.Name = 'PLG custom file plot';
            obj.a = axes;
            obj.a.FontName = 'Arial';
            obj.a.FontSize = 11;
            obj.a.View = [30,30];
            axis vis3d
            obj.a.NextPlot='add';
            obj.a.XLim = [ranger(1,1),ranger(2,1)]; obj.a.XLimMode = 'manual';
            obj.a.YLim = [ranger(1,2),ranger(2,2)]; obj.a.YLimMode = 'manual';
            obj.a.ZLim = [ranger(1,3),ranger(2,3)]; obj.a.ZLimMode = 'manual';
            obj.a.CLim = rangerD; obj.a.CLimMode = 'manual';
            % colour setup
            m = flipud(jet(12)); % reverse jet map
            obj.cMap = colormap(m);
            obj.c = colorbar;
            obj.c.Label.String = 'Diameter';
            obj.c.Label.FontSize = 11;
        end
        function obj = defaultPlot(obj)
            % setup the default data
            
            % vertices - just use a scatter 3
            n = size(obj.vertices,1);
            obj.scVerts = scatter3(...
                obj.vertices(:,1),...
                obj.vertices(:,2),...
                obj.vertices(:,3)...
                );
            % make it pretty
            obj.scVerts.Marker = 'o';
            obj.scVerts.MarkerEdgeColor = 'none';
            obj.scVerts.MarkerFaceColor = 'flat';
            obj.scVerts.SizeData = 7;
            obj.scVerts.CData = obj.sphereDiameter;
            
            % struts - create a graphics array to hold a series of lines
            n = size(obj.struts,1);
            obj.gaStruts = gobjects(n,1);
            nC = size(obj.cMap,1);
            r = obj.a.CLim;
            for inc = 1:n
                points = [obj.vertices(obj.struts(inc,1),:);obj.vertices(obj.struts(inc,2),:)];
                obj.gaStruts(inc) = line(points(:,1),points(:,2),points(:,3));
                % convert strutD to a colour value in cmap
                d = obj.strutDiameter(inc);
                I = round(interp1(r,[1,nC],d));
                obj.gaStruts(inc).Color = obj.cMap(I,:);
            end
        end
    end
end
