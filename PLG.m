classdef PLG
    % The PLG program rewritten in class form
    % only core lattice generation functions will be stored here
    % an input is required to use this
    % EXAMPLE
    % obj = PLG('bcc',12,0.3,...
    %     1,0.3,2,2,2,...
    %     10,10,15,0,0,0);
    % obj = latticeGenerate(obj);
    % save(obj);
    properties (SetAccess=protected)
        latticeType;
        resolution;
        strutDiamter;
        % sphereAddition;
        sphereDiameter;
        vertices;
        faces;
        
        usx;
        usy;
        usz;
        repsx;
        repsy;
        repsz;
        orginx;
        orginy;
        orginz;
        
        validExtensions  = {'xls'; 'csv'; 'custom'};
        filterSpecOut = {'*.stl','3D geometry file';...
            '*.bdf','Nastran input file';...
            '*.inp','Abaqus input file';...
            '*.bin','Binary storage file';...
            '*.xls','Excel format';...
            '*.custom','Custom csv format for debug etc'};
        strutureType;
    end
    methods
        function obj = PLG(varargin)
            % creates the PLG object
            switch numel(varargin)
                case 0
                    % no input no good
                    error('zero inputs is not an option for the PLG class');
                case 1
                    % custom lattice with input as a structure of nodes and
                    % connections
                    obj = loadIn(obj,varargin{1});
                    obj.strutureType = 0;
                case 14
                    %input order same as properties order
                    obj.latticeType = varargin{1};
                    obj.resolution = varargin{2};
                    obj.strutDiamter = varargin{3};
                    if varargin{4}==1
                        obj.sphereDiameter = varargin{5};
                    else
                        obj.sphereDiameter = 0;
                    end
                    
                    obj.usx = varargin{6};
                    obj.usy = varargin{7};
                    obj.usz = varargin{8};
                    obj.repsx = varargin{9};
                    obj.repsy = varargin{10};
                    obj.repsz = varargin{11};
                    obj.orginx = varargin{12};
                    obj.orginy = varargin{13};
                    obj.orginz = varargin{14};
                    
                    obj = latticeGenerate(obj); % generate structure
                    obj = cellReplication(obj); % replicate the unit cells generated above
                    obj = addDiams(obj);
                    obj = cleanLattice(obj); % remove duplicate lattice intersections and struts
                otherwise
                    error('Incorrect number of inputs');
            end
        end
        function obj = latticeGenerate(obj)
            % generates a latticeStructure ready for simulation or saving as an stl etc
            origin = [obj.orginx,obj.orginy,obj.orginz];
            switch obj.latticeType
                case 'bcc' % BCC Cell
                    [obj.vertices, obj.faces] = PLG.bcc(obj.usx,obj.usy,obj.usz,origin);
                case 'bcc2' % BCC Cell with no cantilevers
                    [obj.vertices, obj.faces] = PLG.bcc2(obj.usx,obj.usy,obj.usz,origin);
                case 'bcz' % BCC Cell
                    [obj.vertices, obj.faces] = PLG.bcz(obj.usx,obj.usy,obj.usz,origin);
                case 'fcc'
                    [obj.vertices, obj.faces] = PLG.fcc(obj.usx,obj.usy,obj.usz,origin);
                case 'fccNoXY'
                    [obj.vertices, obj.faces] = PLG.fccNoXY(obj.usx,obj.usy,obj.usz,origin);
                case 'fbcxyz'
                    [obj.vertices, obj.faces] = PLG.fbcxyz(obj.usx,obj.usy,obj.usz,origin);
                case 'fcz'
                    [obj.vertices, obj.faces] = PLG.fcz(obj.usx,obj.usy,obj.usz,origin);
                case 'fbcz'
                    [obj.vertices, obj.faces] = PLG.fbcz(obj.usx,obj.usy,obj.usz,origin);
                case 'bcc_fcc'
                    [obj.vertices, obj.faces] = PLG.bccFcc(obj.usx,obj.usy,obj.usz,origin);
                case 'fbc'
                    [obj.vertices, obj.faces] = PLG.fbc(obj.usx,obj.usy,obj.usz,origin);
                case 'box'
                    [obj.vertices, obj.faces] = PLG.box(obj.usx,obj.usy,obj.usz,origin);
            end
        end
        function obj = cellReplication(obj)
            xPlacement = 0:obj.usx:obj.usx*(obj.repsx-1);
            yPlacement = 0:obj.usy:obj.usy*(obj.repsy-1);
            zPlacement = 0:obj.usz:obj.usz*(obj.repsz-1);
            [XX,YY,ZZ] = ndgrid(xPlacement,yPlacement,zPlacement);
            
            vertOut = arrayfun(@(x,y,z) ...
                [obj.vertices(:,1) + x,...
                obj.vertices(:,2) + y,...
                obj.vertices(:,3) + z]...
                ,XX,YY,ZZ,'UniformOutput',0);
            obj.vertices = cell2mat(vertOut(:));
            
            strutCounter = 0:(numel(XX)-1);
            numStruts = max(max(obj.faces));
            strutCounter = strutCounter*numStruts;
            strutOut = arrayfun(@(x) obj.faces + x,strutCounter,'UniformOutput',0);
            obj.faces = cell2mat(strutOut(:));
        end
        function obj = addDiams(obj)
            % add both sphere and strut diameters
            numNodes = size(obj.vertices,1);
            numStruts = size(obj.faces,1);
            
            obj.sphereDiameter = ones(numNodes,1)*obj.sphereDiameter;
            obj.strutDiamter = ones(numStruts,1)*obj.strutDiamter;
        end
        function obj = loadIn(obj,file)
            parts = strsplit(file,'.');
            extension = parts{end};
            switch extension
                case obj.validExtensions{1}
                    % xls
                    error('TODO');
                otherwise
                    % csv or custom file
                    data = csvread(file);
                    numNodes=data(1,1);
                    numLinks=data(2,1);
                    % data
                    obj.vertices = data(3:numNodes+2,1:3);
                    
                    obj.faces    = data(numNodes+3:numNodes+numLinks+2,1:2);
                    obj.strutDiamter = data(numNodes+3:numNodes+numLinks+2,3);
                    if size(data,2)==3
                        % no sphere diameter supplied
                        obj.sphereDiameter = zeros(numNodes,1);
                    else
                        % sphere diameter supplied
                        obj.sphereDiameter = data(3:numNodes+2,4);
                    end
                    warning('Default resolution of 8 being used');
                    obj.resolution = 8;
            end
        end
        function obj = translate(obj,x,y,z)
            obj.vertices(:,1) = obj.vertices(:,1)+x;
            obj.vertices(:,2) = obj.vertices(:,2)+y;
            obj.vertices(:,3) = obj.vertices(:,3)+z;
        end
        function obj = cleanLattice(obj)
            % cleans the lattice structure including the removal of
            % duplicate vertices and struts
            % remove duplicate faces
            for inc = 1:length(obj.faces)
                ind1 = obj.faces(inc,1);
                ind2 = obj.faces(inc,2);
                if ind1<ind2
                    % faces(inc,1) = ind1;
                    % faces(inc,2) = ind2;
                    % no change
                else
                    obj.faces(inc,1) = ind2;
                    obj.faces(inc,2) = ind1;
                end
            end
            [obj.faces,i] = unique(obj.faces,'rows');
            obj.strutDiamter = obj.strutDiamter(i);
            
            [obj.vertices,i,indexn]=unique(obj.vertices,'rows');
            obj.sphereDiameter = obj.sphereDiameter(i);
            obj.faces = indexn(obj.faces);
        end
        function obj = rotate(obj,wx,wy,wz)
            % rotations are in degrees about the main axes
            thetaX = wx*pi/180;
            thetaY = wy*pi/180;
            thetaZ = wz*pi/180;
            % rotation matricies
            rx = [1           , 0          ,           0;...
                0           , cos(thetaX),-sin(thetaX);...
                0           , sin(thetaX), cos(thetaX)];
            
            ry = [cos(thetaY) ,0           , sin(thetaY);...
                0           ,1           ,0           ;...
                -sin(thetaY),0           , cos(thetaY)];
            
            rz = [cos(thetaZ) ,-sin(thetaZ),           0;...
                sin(thetaZ) , cos(thetaZ),           0;...
                0           ,0           ,           1];
            
            %rotation x then y then zsplit for debugging
            numPoints = length(obj.vertices);
            newPoints = zeros(size(obj.vertices));
            for inc = 1:numPoints
                newPoints(inc,:) = obj.vertices(inc,:)*rx';
            end
            for inc = 1:numPoints
                
            end
            for inc = 1:numPoints
                
            end
            obj.vertices = newPoints;
        end
        function plot(obj)
            % plot the lattice
            f = figure;
            f.Units	= 'Normalized';
            f.Position = [0,0,1,1];
            f.Name = 'STL plotter';
            a = axes;
            a.View = [45,45];
            axis vis3d
            a.NextPlot='add';
            %struts
            p1 = obj.vertices(obj.faces(:,1),:);
            p2 = obj.vertices(obj.faces(:,2),:);
            x = [p1(:,1),p2(:,1)]';
            y = [p1(:,2),p2(:,2)]';
            z = [p1(:,3),p2(:,3)]';
            p = plot3(x,y,z,'Color',[0.3,0.3,0.3,0.5]);
            
            % points
            x = obj.vertices(:,1);
            y = obj.vertices(:,2);
            z = obj.vertices(:,3);
            s = scatter3(x,y,z,'MarkerFaceColor',[0.9,0.5,0],'MarkerEdgeColor','none');
            
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
    end
    methods % save out methods
        function save(obj)
            % this overloads the save function and allows saving out in various
            % formats however this only saves the latticeStructure structure
            [fileName,pathName,filterIndex] = uiputfile(obj.filterSpecOut);
            switch filterIndex
                case 1
                    % stl
                    saveStl(obj,fileName,pathName);
                case 2
                    % Nastran
                    error('TODO');
                case 3
                    % ABAQUS
                    error('TODO');
                case 4
                    % BINARY
                    error('TODO');
                case 5
                    saveExcel(obj,fileName,pathName)
                case 6
                    saveCustom(obj,fileName,pathName)
                otherwise
                    error('No file saved')
            end
        end
        function saveStl(obj,fileName,pathName)
            fullName  = [pathName,fileName];
            numFacets = obj.resolution*4;%number of facets created for one strut
            numLinks = size(obj.faces,1);
            numVertices = size(obj.vertices,1);
            
            totalFacetsNoBall = numFacets*numLinks;
            totalFacetsWithBall = totalFacetsNoBall + 8^2*(8-1)/4*numVertices;
            
            fid=fopen(fullName,'w');
            fprintf(fid, '%-80s', 'fast stl generator'); %binary write information
            
            if ~isempty(obj.sphereDiameter)
                fwrite(fid,uint32(totalFacetsWithBall),'uint32'); %stl binary header file contains the total number of facets in the stl file
                fid = ballCreate(obj,fid);
            else
                fwrite(fid,uint32(totalFacetsNoBall),'uint32'); %stl binary header file contains the total number of facets in the stl file
            end
            fid = faceCreate(obj,fid);
            fclose(fid);
        end
        function saveExcel(obj,fileName,pathName)
            % saves data to an excel compatible format
            fullName = [pathName,fileName];
            
            numNodes=size(obj.vertices,1);
            numLinks=size(obj.faces,1);
            xlswrite(fullName,numNodes,'Sheet1','A1');
            xlswrite(fullName,numLinks,'Sheet1','A2');
            
            locationVertices=strcat('A3:C',num2str(numNodes+2));
            xlswrite(fullName,obj.vertices,'Sheet1',locationVertices);
            locationFaces=strcat('A',num2str(numNodes+3),':B',num2str(numNodes+2+numLinks));
            xlswrite(fullName,obj.faces,'Sheet1',locationFaces);
            
            % optional write depending on load type
            if obj.strutureType ==0
                % loaded custom
                locationSpheres = strcat('D3:D',num2str(numNodes+2));
                xlswrite(fullName,obj.spheres,'Sheet1',locationSpheres);
                locationDiameters = strcat('C',num2str(numNodes+3),':C',num2str(numNodes+2+numLinks));
                xlswrite(fullName,obj.diameters,'Sheet1',locationDiameters);
            else
                % loaded regular
                if ~isempty(obj.sphereDiameter)
                    locationSpheres = strcat('D3:D',num2str(numNodes+2));
                    data = ones(numNodes,1)*obj.sphereDiameter;
                    xlswrite(fullName,data,'Sheet1',locationSpheres);
                else
                    locationSpheres = strcat('D3:D',num2str(numNodes+2));
                    data = ones(numNodes,1)*0;
                    xlswrite(fullName,data,'Sheet1',locationSpheres);
                end
                locationDiameters = strcat('C',num2str(numNodes+3),':C',num2str(numNodes+2+numLinks));
                data = ones(numLinks,1)*obj.strutDiamter;
                xlswrite(fullName,data,'Sheet1',locationDiameters);
            end
        end
        function saveCustom(obj,fileName,pathName)
            % saves data to an excel compatible format
            fullName = [pathName,fileName];
            
            numNodes=size(obj.vertices,1);
            numLinks=size(obj.faces,1);
            
            dlmwrite(fullName,numNodes);
            dlmwrite(fullName,numLinks,'-append');
            data = [obj.vertices,obj.sphereDiameter];
            dlmwrite(fullName,data,'-append');
            
            data = [obj.faces,obj.strutDiamter];
            dlmwrite(fullName,data,'-append');
        end
        function fid = faceCreate(obj,fid)
            radius = obj.strutDiamter/2;
            numLinks = length(obj.faces);
            
            % calculate points on each facet then write said points to file
            % in the correct order
            for i=1:numLinks
                point1 = obj.vertices(obj.faces(i,1),:);
                point2 = obj.vertices(obj.faces(i,2),:);
                vector = point2-point1;
                u1 = vector/norm(vector);
                if u1(3)==1
                    v1 = [1,0,0];
                else
                    v1 = cross([0,0,1],u1);
                    v1 = v1/norm(v1);
                end
                offset = radius(i)*v1;
                
                vert1end = zeros(obj.resolution,3);
                vert2end = zeros(obj.resolution,3);
                for j=1:obj.resolution
                    Qrot1 = PLG.qGetRotQuaternion((j-1)*2*pi/obj.resolution, u1);
                    absolutePointRotation = PLG.qRotatePoint(offset', Qrot1)';
                    % end 1
                    vert1end(j,:)=absolutePointRotation+point1;
                    % end 2
                    vert2end(j,:)=absolutePointRotation+point2;
                end
                % scatter3(vertOut(:,1),vertOut(:,2),vertOut(:,3)) % scatter
                %% join faces
                for j=1:obj.resolution
                    % end of strut at point 1
                    datOut = circshift(vert1end,j);
                    facet_a=point1;
                    facet_b=datOut(1,:);
                    facet_c=datOut(2,:);
                    normal=cross(facet_b-facet_a,facet_a-facet_a);
                    fwrite(fid,normal,'float32');           % write normal vector floating point numbers
                    fwrite(fid,facet_c,'float32');   % first vertex (x,y,z)
                    fwrite(fid,facet_b,'float32');   % second vertex
                    fwrite(fid,facet_a,'float32');   % Third vertex
                    fwrite(fid,0,'uint16','l');
                    % end of strut at point 2
                    datOut = circshift(vert2end,j);
                    facet_a=point2;
                    facet_b=datOut(1,:);
                    facet_c=datOut(2,:);
                    normal=cross(facet_b-facet_a,facet_a-facet_a);
                    fwrite(fid,normal,'float32');           % write normal vector floating point numbers
                    fwrite(fid,facet_c,'float32');   % first vertex (x,y,z)
                    fwrite(fid,facet_b,'float32');   % second vertex
                    fwrite(fid,facet_a,'float32');   % Third vertex
                    fwrite(fid,0,'uint16','l');
                    % along direction point 1 to point 2
                    datOut1 = circshift(vert1end,j);
                    datOut2 = circshift(vert2end,j);
                    facet_a=datOut1(1,:);
                    facet_b=datOut2(1,:);
                    facet_c=datOut2(2,:);
                    normal=cross(facet_b-facet_a,facet_a-facet_a);
                    fwrite(fid,normal,'float32');           % write normal vector floating point numbers
                    fwrite(fid,facet_c,'float32');   % first vertex (x,y,z)
                    fwrite(fid,facet_b,'float32');   % second vertex
                    fwrite(fid,facet_a,'float32');   % Third vertex
                    fwrite(fid,0,'uint16','l');
                    % along direction point 2 to point 1
                    datOut1 = circshift(vert1end,j);
                    datOut2 = circshift(vert2end,j);
                    facet_a=datOut1(1,:);
                    facet_b=datOut1(2,:);
                    facet_c=datOut2(2,:);
                    normal=cross(facet_b-facet_a,facet_a-facet_a);
                    fwrite(fid,normal,'float32');           % write normal vector floating point numbers
                    fwrite(fid,facet_a,'float32');   % first vertex (x,y,z)
                    fwrite(fid,facet_b,'float32');   % second vertex
                    fwrite(fid,facet_c,'float32');   % Third vertex
                    fwrite(fid,0,'uint16','l');
                end
            end
        end
        function fid = ballCreate(obj,fid)
                [x,y,z]=sphere(8); %create sphere with higher accuracy
                ball.faces= convhull([x(:), y(:), z(:)]); %create triangle links
                sizer = size(ball.faces,1);
                ball.vertices=[x(:),y(:),z(:)]; %store the points
                
                for i=1:size(obj.vertices,1)
                    offset=ball.vertices*obj.sphereDiameter(i)/2;
                    target=[offset(:,1)+obj.vertices(i,1),offset(:,2)+obj.vertices(i,2),offset(:,3)+obj.vertices(i,3)];
                    for j=1:sizer
                        %get values first end
                        facet_a=target(ball.faces(j,1),:);
                        facet_b=target(ball.faces(j,2),:);
                        facet_c=target(ball.faces(j,3),:);
                        normal=cross(facet_b-facet_a,facet_c-facet_a);
                        %write values
                        fwrite(fid,normal,'float32');           % write normal vector floating point numbers
                        fwrite(fid,facet_a,'float32');   % first vertex (x,y,z)
                        fwrite(fid,facet_b,'float32');   % second vertex
                        fwrite(fid,facet_c,'float32');   % Third vertex
                        fwrite(fid,32767,'uint16','l');
                    end
                end
        end
    end
    methods (Static) %cell type methods
        function [nodes, faces]= bcc(sx, sy, sz, origin)
            nodes=[origin;(origin(1)+sx),origin(2),origin(3);origin(1)+sx,origin(2)+sy,origin(3);...
                origin(1),origin(2)+sy,origin(3);origin(1),origin(2),origin(3)+sz;...
                origin(1)+sx,origin(2),origin(3)+sz;origin(1)+sx,origin(2)+sy,origin(3)+sz;...
                origin(1),origin(2)+sy,origin(3)+sz;...
                origin(1)+sx/2,origin(2)+sy/2,origin(3)+sz/2];
            faces=[1,9;2,9;3,9;4,9;5,9;6,9;7,9;8,9];
        end
        function [nodes, faces]= bcz(sx, sy, sz, origin)
            nodes=[origin;...                                        1
                origin(1)+sx,origin(2),origin(3);...              2
                origin(1)+sx,origin(2)+sy,origin(3);...           3
                origin(1),origin(2)+sy,origin(3);...              4
                origin(1),origin(2),origin(3)+sz;...              5
                origin(1)+sx,origin(2),origin(3)+sz;...           6
                origin(1)+sx,origin(2)+sy,origin(3)+sz;...        7
                origin(1),origin(2)+sy,origin(3)+sz;...           8
                origin(1)+sx/2,origin(2)+sy/2,origin(3)+sz/2];  % 9
            faces=[1,9;...
                1,5;...
                2,9;...
                2,6;...
                3,9;...
                3,7;...
                4,9;...
                4,8;...
                5,9;...
                6,9;...
                7,9;....
                8,9];
        end
        function [nodes, faces]= bcc2(sx, sy, sz, origin)
            nodes=[origin(1),origin(2)+sy/2,origin(3)+sz/2;...
                origin(1)+sx,origin(2)+sy/2,origin(3)+sz/2;...
                origin(1)+sx/2,origin(2),origin(3);...
                origin(1)+sx/2,origin(2),origin(3)+sz;...
                origin(1)+sx/2,origin(2)+sy,origin(3)+sz;...
                origin(1)+sx/2,origin(2)+sy,origin(3)];
            faces=[3,1;...
                4,1;...
                5,1;...
                6,1;...
                2,3;...
                2,4;...
                2,5;...
                2,6];
        end
        function [nodes, faces]= fcc(sx, sy, sz, origin)
            nodes=[origin;...                                           1  0
                origin(1)+sx,origin(2),origin(3);...                 2  0
                origin(1)+sx,origin(2)+sy,origin(3);...              3  0
                origin(1),origin(2)+sy,origin(3);...                 4  0
                origin(1),origin(2),origin(3)+sz;...                 5  1
                origin(1)+sx,origin(2),origin(3)+sz;...              6  1
                origin(1)+sx,origin(2)+sy,origin(3)+sz;...           7  1
                origin(1),origin(2)+sy,origin(3)+sz;...              8  1
                origin(1)+sx/2,origin(2)+sy/2,origin(3);...          9  0
                origin(1)+sx/2,origin(2)+sy/2,origin(3)+sz;...       10 1
                origin(1)+sx/2,origin(2),origin(3)+sz/2;...          11 2
                origin(1)+sx/2,origin(2)+sy,origin(3)+sz/2;...       12 2
                origin(1),origin(2)+sy/2,origin(3)+sz/2;...          13 2
                origin(1)+sx,origin(2)+sy/2,origin(3)+sz/2]; %       14 2
            faces=[1,13;...
                1,9;...
                1,11;...
                2,9;...
                2,14;...
                2,11;...
                3,12;...
                3,14;...
                3,9;...
                4,12;...
                4,9;...
                4,13;...
                5,13;...
                5,11;...
                5,10;...
                6,10;...
                6,14;...
                6,11;...
                7,10;...
                7,14;...
                7,12;...
                8,12;...
                8,13;...
                8,10];
        end
        function [nodes, faces]= fccNoXY(sx, sy, sz, origin)
            nodes=[origin;...                                        1  0
                origin(1)+sx,origin(2),origin(3);...                 2  0
                origin(1)+sx,origin(2)+sy,origin(3);...              3  0
                origin(1),origin(2)+sy,origin(3);...                 4  0
                origin(1),origin(2),origin(3)+sz;...                 5  1
                origin(1)+sx,origin(2),origin(3)+sz;...              6  1
                origin(1)+sx,origin(2)+sy,origin(3)+sz;...           7  1
                origin(1),origin(2)+sy,origin(3)+sz;...              8  1
                origin(1)+sx/2,origin(2),origin(3)+sz/2;...          9  2
                origin(1)+sx/2,origin(2)+sy,origin(3)+sz/2;...       10 2
                origin(1),origin(2)+sy/2,origin(3)+sz/2;...          11 2
                origin(1)+sx,origin(2)+sy/2,origin(3)+sz/2]; %       12 2
            faces=[1,11;...
                1,9;...
                2,12;...
                2,9;...
                3,10;...
                3,12;...
                4,10;...
                4,11;...
                5,11;...
                5,9;...
                6,12;...
                6,9;...
                7,12;...
                7,10;...
                8,10;...
                8,11];
        end
        function [nodes, faces]= fcz(sx, sy, sz, origin)
            nodes=[origin;...
                origin(1)+sx,origin(2),origin(3);...
                origin(1)+sx,origin(2)+sy,origin(3);...
                origin(1),origin(2)+sy,origin(3);...
                origin(1),origin(2),origin(3)+sz;...
                origin(1)+sx,origin(2),origin(3)+sz;...
                origin(1)+sx,origin(2)+sy,origin(3)+sz;...
                origin(1),origin(2)+sy,origin(3)+sz;...
                origin(1)+sx/2,origin(2),origin(3)+sz/2;...
                origin(1)+sx,origin(2)+sy/2,origin(3)+sz/2;...
                origin(1)+sx/2,origin(2)+sy,origin(3)+sz/2;...
                origin(1),origin(2)+sy/2,origin(3)+sz/2];
            faces=[1,5;1,9;1,12;2,9;2,6;2,10;3,10;3,11;3,7;4,11;4,12;4,8;5,9;5,12;6,9;6,10;7,11;7,10;8,11;8,12];
        end
        function [nodes, faces]= box(sx, sy, sz, origin)
            nodes=[origin;(origin(1)+sx),origin(2),origin(3);origin(1)+sx,origin(2)+sy,origin(3);...
                origin(1),origin(2)+sy,origin(3);origin(1),origin(2),origin(3)+sz;...
                origin(1)+sx,origin(2),origin(3)+sz;origin(1)+sx,origin(2)+sy,origin(3)+sz;...
                origin(1),origin(2)+sy,origin(3)+sz];
        end
        function [nodes, faces]= bccFcc(sx, sy, sz, origin)
            nodes=[origin;(origin(1)+sx),origin(2),origin(3);...
                origin(1)+sx,origin(2)+sy,origin(3);...
                origin(1),origin(2)+sy,origin(3);origin(1),origin(2),origin(3)+sz;...
                origin(1)+sx,origin(2),origin(3)+sz;...
                origin(1)+sx,origin(2)+sy,origin(3)+sz;...
                origin(1),origin(2)+sy,origin(3)+sz;...
                origin(1)+sx/2,origin(2)+sy/2,origin(3)+sz/2;...
                origin(1)+sx/2,origin(2)+sy/2,origin(3);...
                origin(1)+sx/2,origin(2)+sy/2,origin(3)+sz;...
                origin(1)+sx/2,origin(2),origin(3)+sz/2;...
                origin(1)+sx/2,origin(2)+sy,origin(3)+sz/2;...
                origin(1),origin(2)+sy/2,origin(3)+sz/2;...
                origin(1)+sx,origin(2)+sy/2,origin(3)+sz/2];
        end
        function [nodes, faces]= fbc(sx, sy, sz, origin)
            nodes=[origin;...
                origin(1)+sx,origin(2),origin(3);...
                origin(1),origin(2)+sy,origin(3);...
                origin(1)+sx,origin(2)+sy,origin(3);...
                origin(1),origin(2),origin(3)+sz;...
                origin(1)+sx,origin(2),origin(3)+sz;...
                origin(1),origin(2)+sy,origin(3)+sz;...
                origin(1)+sx,origin(2)+sy,origin(3)+sz];
            faces=[2,8;4,6;1,7;3,5;1,6;2,5;4,7;3,8;... faces
                1,8;2,7;4,5;3,6]; %body
        end
        function [nodes, faces]= fbcxyz(sx, sy, sz, origin)
            nodes=[origin;...
                origin(1)+sx,origin(2),origin(3);...
                origin(1),origin(2)+sy,origin(3);...
                origin(1)+sx,origin(2)+sy,origin(3);...
                origin(1),origin(2),origin(3)+sz;...
                origin(1)+sx,origin(2),origin(3)+sz;...
                origin(1),origin(2)+sy,origin(3)+sz;...
                origin(1)+sx,origin(2)+sy,origin(3)+sz];
            faces=[1,2;1,3;4,3;4,2;2,6;6,8;8,4;1,5;5,7;7,3;5,6;7,8;... edges
                1,4;2,3;2,8;4,6;1,7;3,5;1,6;2,5;4,7;3,8;7,6;5,8;... faces
                1,8;2,7;4,5;3,6]; %body
        end
        function [nodes, faces]= fbcz(sx, sy, sz, origin)
            nodes=[origin;...
                origin(1)+sx,origin(2),origin(3);...
                origin(1),origin(2)+sy,origin(3);...
                origin(1)+sx,origin(2)+sy,origin(3);...
                origin(1),origin(2),origin(3)+sz;...
                origin(1)+sx,origin(2),origin(3)+sz;...
                origin(1),origin(2)+sy,origin(3)+sz;...
                origin(1)+sx,origin(2)+sy,origin(3)+sz];
            faces=[2,6;8,4;1,5;7,3;... edges
                2,8;4,6;1,7;3,5;1,6;2,5;4,7;3,8;... faces
                1,8;2,7;4,5;3,6]; %body
        end
    end
    methods (Static) % quaternian methods
        function Q = qConj( Q1 )
            % qConj: quaternion conjugation
            Q1 = reshape( Q1, 4, 1 );
            Q = [Q1(1);-Q1(2:4)];
        end
        function Q = qGetRotQuaternion( teta, rotationVector )
            % qCreateRotQuaternion: outputs a quaternion which is used to perform
            % rotation
            norm = sqrt(sum( rotationVector .* rotationVector ));
            if( numel( rotationVector ) ~= 3 )
                fprintf( 'rotationVector should have 3 coordinates!\r' );
                return;
            end
            
            rotationVector = reshape( rotationVector, 3, 1 );
            if( norm > 0 )
                v = rotationVector / norm;
                Q = [ cos( teta/2 ); v*sin( teta/2 )];
            else
                fprintf( 'rotationVector cannot be 0\n\r' );
            end
        end
        function Q = qInv( Q1 )
            % qInv: quaternion reciprocal (inverse)
            % Q = qInv( Q1 )
            Q = PLG.qConj( Q1 ) ./ PLG.qLength( Q1 )^2;
        end
        function d = qLength( Q )
            % qLength: quaternion length (norm)
            d = sqrt( sum( Q .* Q ));
        end
        function Q = qMul( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10 )
            % qMul: quaternion multiplication
            if( nargin >= 2 )
                Q = PLG.qMul2( Q1, Q2 );
            end
            if( nargin >= 3 )
                Q = PLG.qMul2( Q, Q3 );
            end
            if( nargin >= 4 )
                Q = PLG.qMul2( Q, Q4 );
            end
            if( nargin >= 5 )
                Q = PLG.qMul2( Q, Q5 );
            end
            if( nargin >= 6 )
                Q = PLG.qMul2( Q, Q6 );
            end
            if( nargin >= 7 )
                Q = PLG.qMul2( Q, Q7 );
            end
            if( nargin >= 8 )
                Q = PLG.qMul2( Q, Q8 );
            end
            if( nargin >= 9 )
                Q = PLG.qMul2( Q, Q9 );
            end
            if( nargin >= 10)
                Q = PLG.qMul2( Q, Q10);
            end
        end
        function Q = qMul2( Q1, Q2 )
            % qMul: quaternion multiplication
            
            s1 = Q1(1);
            s2 = Q2(1);
            v1 = Q1(2:4);
            v2 = Q2(2:4);
            
            s =s1*s2 - dot( v1,v2);
            v = s1*v2 + s2*v1 + cross( v1, v2 );
            v = reshape( v, 3, 1 );
            Q = [s;v];
        end
        function Protated = qRotatePoint( P, Qrotation )
            % qRotatePoint: rotate a point according to rotation quaternion
            Q1 = [ 0; P ];
            Q = PLG.qMul( Qrotation, Q1, PLG.qInv( Qrotation ) );
            Protated = Q(2:4);
        end
    end
end
