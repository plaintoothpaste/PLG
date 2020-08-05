classdef plotPLG
    % The plotPLG program will plot custom files with desired settings
    properties (SetAccess=protected)
        f
        a
        c
        colour_map % colourMap
        h_struts % a grapics array for the struts
        h_verts % a scatter3 of the verts
    end
    properties (SetAccess=private)
        vertices
        struts
        strut_dia
        sphere_dia
    end
    methods
        function obj = plotPLG(fName)
            % load the data and prep a basic plot
            obj = load(obj,fName);
            obj = figureSetupAndPosition(obj);
            obj = defaultPlot(obj);
        end
        function obj = setColor2incline(obj)
            % sets the color map to a default inclination option
        end
        function obj = setColor2diameter(obj)
            % sets color map to a default diameter option
            dia_range = [min(d), max(d)];
            if dia_range(1)==dia_range(2)% in case value 1 and 2 match
                dia_range(1) = dia_range(1)-0.1; 
                dia_range(2) = dia_range(2)-0.1;
            end
            
            obj.a.CLim = dia_range; obj.a.CLimMode = 'manual';
            % colour setup
            m = flipud(jet(12)); % reverse jet map
            obj.colour_map = colormap(m);
            obj.c.Label.String = 'Diameter';
        end
    end
    methods (Access=private) % called when plotPLG is used
        function obj = load(obj,file)
            % load a custom beam input file to generate a lattice structure
            [~,~,ext] = fileparts(file);
            if ~strcmp(ext,'.custom')
                error('Load in format must be .custom');
            end
            
            data = csvread(file);
            num_nodes=data(1,1);
            num_links=data(2,1);
            vert_loc = 3:num_nodes+2;
            strut_loc = num_nodes+3:num_nodes+num_links+2;
            
            obj.vertices = data(vert_loc,1:3);
            obj.struts   = data(strut_loc,1:2);
            obj.strut_dia = data(strut_loc,3);
            if size(data,2)==3
                obj.sphere_dia = zeros(num_nodes,1);% no sphere diameter supplied
            else
                obj.sphere_dia = data(vert_loc,4);
            end
        end
        function obj = figureSetupAndPosition(obj)
            % create the figure (f) and axes (a) and colour bar
            vertices_range = [min(obj.vertices);max(obj.vertices)];
            
            % create a standard figure 3:2 ratio
            obj.f = figure;
            obj.f.Units	= 'centimeters';
            obj.f.Position = [3,3,18,12];
            obj.f.Name = 'PLG custom file plot';
            
            % create an axes with the limits== extents of the lattice
            obj.a = axes;
            obj.a.FontName = 'Arial';
            obj.a.FontSize = 11;
            obj.a.View = [30,30];
            axis vis3d
            obj.a.NextPlot='add';
            obj.a.XLim = [vertices_range(1,1),vertices_range(2,1)]; 
            obj.a.YLim = [vertices_range(1,2),vertices_range(2,2)]; 
            obj.a.ZLim = [vertices_range(1,3),vertices_range(2,3)];
            obj.a.XLimMode = 'manual'; obj.a.YLimMode = 'manual'; obj.a.ZLimMode = 'manual';
            % pre setup of color bar
            obj.c = colorbar;
            obj.c.Label.FontSize = 11;
            
            % use diamter colour plot as the default.
            obj = setColor2diameter(obj);
        end
        function obj = defaultPlot(obj)
            % setup the default data
            
            % vertices - just use a scatter 3
            n = size(obj.vertices,1);
            obj.h_verts = scatter3(...
                obj.vertices(:,1),...
                obj.vertices(:,2),...
                obj.vertices(:,3)...
                );
            % make it pretty
            obj.h_verts.Marker = 'o';
            obj.h_verts.MarkerEdgeColor = 'none';
            obj.h_verts.MarkerFaceColor = 'flat';
            obj.h_verts.SizeData = 7;
            obj.h_verts.CData = obj.sphere_dia;
            
            % struts - create a graphics array to hold a series of lines
            n = size(obj.struts,1);
            obj.h_struts = gobjects(n,1);
            nC = size(obj.colour_map,1);
            r = obj.a.CLim;
            for inc = 1:n
                points = [obj.vertices(obj.struts(inc,1),:);obj.vertices(obj.struts(inc,2),:)];
                obj.h_struts(inc) = line(points(:,1),points(:,2),points(:,3));
                % convert strutD to a colour value in colour_map
                d = obj.strut_dia(inc);
                I = round(interp1(r,[1,nC],d));
                obj.h_struts(inc).Color = obj.colour_map(I,:);
            end
        end
    end
end
