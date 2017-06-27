classdef PLG
    % The PLG program rewritten in class form
    % only core lattice generation functions will be stored here
    % an input is required to use this
    % EXAMPLE
    % obj = PLG('bcc',12,0.3,...
    %     1,0.5,2,2,2,...
    %     4,4,5,0,0,0);
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
        tolerance; % defined as 1/1000 of the shortest length present
        dx; % length of a strut
    end
    methods
        function obj = PLG(varargin)
            % creates the PLG object
            switch numel(varargin)
                case 0
                    % no input no good
                    error('zero inputs is not an option for the PLG class');
                case 1
                    % import a custom lattice file containing beam and node
                    % definitions see load function for more information
                    obj = load(obj,varargin{1});
                    obj.strutureType = 0;
                case 14
                    %Generate a new regular lattice
                    %   input order same as properties order
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
                    obj = addDiams(obj);        % apply unique diameters to each strut
                    obj = cleanLattice(obj);    % remove duplicate lattice intersections and struts
                otherwise
                    error('Incorrect number of inputs');
            end
        end
        function obj = latticeGenerate(obj)
            % generates a latticeStructure ready for simulation or saving as an stl etc
            % all the below unit cells are stored in their own method section
            switch obj.latticeType
                case 'bcc' % BCC Cell
                    [obj.vertices, obj.faces] = PLG.bcc(obj.usx,obj.usy,obj.usz);
                case 'bcc2' % BCC Cell with no cantilevers
                    [obj.vertices, obj.faces] = PLG.bcc2(obj.usx,obj.usy,obj.usz);
                case 'bcz' % BCC Cell
                    [obj.vertices, obj.faces] = PLG.bcz(obj.usx,obj.usy,obj.usz);
                case 'fcc'
                    [obj.vertices, obj.faces] = PLG.fcc(obj.usx,obj.usy,obj.usz);
                case 'fccNoXY'
                    [obj.vertices, obj.faces] = PLG.fccNoXY(obj.usx,obj.usy,obj.usz);
                case 'fbcxyz'
                    [obj.vertices, obj.faces] = PLG.fbcxyz(obj.usx,obj.usy,obj.usz);
                case 'fcz'
                    [obj.vertices, obj.faces] = PLG.fcz(obj.usx,obj.usy,obj.usz);
                case 'fbcz'
                    [obj.vertices, obj.faces] = PLG.fbcz(obj.usx,obj.usy,obj.usz);
                case 'bcc_fcc'
                    [obj.vertices, obj.faces] = PLG.bccFcc(obj.usx,obj.usy,obj.usz);
                case 'fbc'
                    [obj.vertices, obj.faces] = PLG.fbc(obj.usx,obj.usy,obj.usz);
                case 'box'
                    [obj.vertices, obj.faces] = PLG.box(obj.usx,obj.usy,obj.usz);
            end
            % translsate the unit cell to the correct location
            obj = translate(obj,obj.orginx,obj.orginy,obj.orginz);
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
        function obj = load(obj,file)
            % load a custom beam input file to generate a lattice structure
            parts = strsplit(file,'.');
            extension = parts{end};
            switch extension
                case obj.validExtensions{1}
                    % xls
                    error('TODO');
                case obj.validExtensions{2}
                    % csv
                    data = csvread(file);
                    numNodes=data(1,1);
                    numLinks=data(2,1);
                    % data
                    obj.vertices = data(3:numNodes+2,1:3);
                    
                    obj.faces    = data(numNodes+3:numNodes+numLinks+2,1:2);
                    obj.strutDiameter = data(numNodes+3:numNodes+numLinks+2,3);
                    if size(data,2)==3
                        % no sphere diameter supplied
                        obj.sphereDiameter = zeros(numNodes,1);
                    else
                        % sphere diameter supplied
                        obj.sphereDiameter = data(3:numNodes+2,4);
                    end
                case obj.validExtensions{3}
                    % custom
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
            end
        end
        function obj = translate(obj,x,y,z)
            obj.vertices(:,1) = obj.vertices(:,1)+x;
            obj.vertices(:,2) = obj.vertices(:,2)+y;
            obj.vertices(:,3) = obj.vertices(:,3)+z;
        end
        function obj = cleanLattice(obj)
            % cleans the lattice structure including the removal of
            % duplicate vertices and removes duplicate faces
            
            % duplicate faces
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
            
            % duplicate vertices
            [obj.vertices,i,indexn]=uniquetol(obj.vertices,obj.tolerance,'ByRows',1,'DataScale',1);
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
            % plot the lattice with nodes highlighted
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
        function obj = splitStruts(obj)
            % takes a badly defined input lattice and scoures it for any struts that intersect and
            % but do not connect to each other. It then breaks these struts and creates a new clean
            % lattice
            % this functions should only be used if the user knows what they are doing and may
            % destroy unique diameter data
            maxInd = 0; % current highest index in facesOut
            lengthFaces = 0; % current length of facesOut
            numFaces = length(obj.faces);
            vertsOut = [];
            facesOut = [];
            splitStruts = cell(numFaces,1);
            % get tol
            obj = getTolerance(obj);
            tol = obj.tolerance;
            struts2Check = 1:numFaces;
            boundingBox = PLG.getBound(obj.vertices,obj.faces,struts2Check);
            %% main loop
            for inc=1:numFaces
                % all struts with intersecting BB
                struts2CheckTmp = struts2Check;
                struts2CheckTmp(inc) = [];
                currentBB = boundingBox{inc};
                currentP1 = obj.vertices(obj.faces(inc,1),:);
                currentP2 = obj.vertices(obj.faces(inc,2),:);
                isIntersect = PLG.findBbIntersect(currentBB,boundingBox(struts2CheckTmp));
                potentialStruts = struts2CheckTmp(isIntersect);
                % remove that are already connected properly
                potentialStruts = PLG.removeNormalConStruts(currentP1,currentP2,potentialStruts,obj.faces,obj.vertices,tol);
                potentialStruts = PLG.removeParralelStruts(currentP1,currentP2,potentialStruts,obj.faces,obj.vertices);
                % check if any faces are already split if so grab their sub BB and
                % eleminate any that do not intersect
                % potentialStruts = PLG.findSplitStrutBB(splitStruts,potentialStruts,currentBB,currentP1,currentP2,facesOut,vertsOut,tol);
                potentialStruts = [potentialStruts',ones(length(potentialStruts),1)];
                if isempty(potentialStruts)
                    % no struts to break up
                    vertsOut = [vertsOut;currentP1;currentP2];
                    facesOut = [facesOut;maxInd+1,maxInd+2];
                    maxInd = maxInd+2;
                    lengthFaces = lengthFaces+1; % number of new faces
                else
                    [vertsTmpOut,facesTmpOut] = PLG.findIntersect(currentP1,currentP2,...
                        potentialStruts,obj.faces,obj.vertices,...
                        facesOut,vertsOut,tol);
                    vertsOut = [vertsOut;vertsTmpOut];
                    facesTmpOut = facesTmpOut+maxInd;
                    % debug: 
                    % patch('Faces',obj.faces(inc,:),'Vertices',obj.vertices,'EdgeColor',[1,0,0]);
                    % patch('Faces',obj.faces(potentialStruts,:),'Vertices',obj.vertices); hold on;
                    % scatter3(vertsOut(:,1),vertsOut(:,2),vertsOut(:,3));

                    maxInd = facesTmpOut(end,2);
                    facesOut = [facesOut;facesTmpOut];
                    st = lengthFaces+1; % number of new faces
                    en = st+size(facesTmpOut,1)-1;
                    splitStruts{inc} = [st,en];
                    lengthFaces = en;
                end
                
            end
            %% final assignment and cleanup note that this destroys strut diameter
            obj.vertices = vertsOut;
            obj.faces = facesOut;
            % remake sphere and strut diams same as the first value
            obj.strutDiamter = obj.strutDiamter(1);
            obj.sphereDiameter = obj.sphereDiameter(1);
            obj = addDiams(obj);
            obj = cleanLattice(obj);
        end
        function obj = getTolerance(obj)
            % get the shortest strut and divide by 1000
            verts1 = obj.vertices(obj.faces(:,1),:);
            verts2 = obj.vertices(obj.faces(:,2),:);
            lengthVerts = sum(sqrt((verts1-verts2).^2),2);
            obj.tolerance = min(lengthVerts)/10;
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
                    % TODO
                case 3
                    % ABAQUS
                    % TODO
                case 4
                    % BINARY
                    % TODO
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
        function [nodes, faces]= bcc(sx, sy, sz)
            nodes=[0 , 0, 0;...
                sx, 0, 0;...
                sx,sy, 0;...
                0 ,sy, 0;...
                0 , 0,sz;...
                sx, 0,sz;...
                sx,sy,sz;...
                0 ,sy,sz;...
                sx/2,sy/2,sz/2];
            faces=[1,9;2,9;3,9;4,9;5,9;6,9;7,9;8,9];
        end
        function [nodes, faces]= bcz(sx, sy, sz)
            nodes=[0 , 0, 0;...  1
                sx, 0, 0;...  2
                sx,sy, 0;...  3
                0 ,sy, 0;...  4
                0 , 0,sz;...  5
                sx,0,sz;...           6
                sx,sy,sz;...        7
                0,sy,sz;...           8
                sx/2,sy/2,sz/2];  % 9
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
        function [nodes, faces]= bcc2(sx, sy, sz)
            nodes=[0,sy/2,sz/2;...
                sx,sy/2,sz/2;...
                sx/2,0,0;...
                sx/2,0,sz;...
                sx/2,sy,sz;...
                sx/2,sy,0];
            faces=[3,1;...
                4,1;...
                5,1;...
                6,1;...
                2,3;...
                2,4;...
                2,5;...
                2,6];
        end
        function [nodes, faces]= fcc(sx, sy, sz)
            nodes=[0,0,0;...               1  0
                sx,0,0;...                 2  0
                sx,sy,0;...              3  0
                0,sy,0;...                 4  0
                0,0,sz;...                 5  1
                sx,0,sz;...              6  1
                sx,sy,sz;...           7  1
                0,sy,sz;...              8  1
                sx/2,sy/2,0;...          9  0
                sx/2,sy/2,sz;...       10 1
                sx/2,0,sz/2;...          11 2
                sx/2,sy,sz/2;...       12 2
                0,sy/2,sz/2;...          13 2
                sx,sy/2,sz/2]; %       14 2
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
        function [nodes, faces]= fccNoXY(sx, sy, sz)
            nodes=[0,0,0;...             1  0
                sx,0,0;...                 2  0
                sx,sy,0;...              3  0
                0,sy,0;...                 4  0
                0,0,sz;...                 5  1
                sx,0,sz;...              6  1
                sx,sy,sz;...           7  1
                0,sy,sz;...              8  1
                sx/2,0,sz/2;...          9  2
                sx/2,sy,sz/2;...       10 2
                0,sy/2,sz/2;...          11 2
                sx,sy/2,sz/2]; %       12 2
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
        function [nodes, faces]= fcz(sx, sy, sz)
            nodes=[0,0,0;...
                sx,0,0;...
                sx,sy,0;...
                0,sy,0;...
                0,0,sz;...
                sx,0,sz;...
                sx,sy,sz;...
                0,sy,sz;...
                sx/2,0,sz/2;...
                sx,sy/2,sz/2;...
                sx/2,sy,sz/2;...
                0,sy/2,sz/2];
            faces=[1,5;1,9;1,12;2,9;2,6;2,10;3,10;3,11;3,7;4,11;4,12;4,8;5,9;5,12;6,9;6,10;7,11;7,10;8,11;8,12];
        end
        function [nodes, faces]= box(sx, sy, sz)
            nodes=[0,0,0;... 1
                sx,0,0;... 2
                sx,sy,0;... 3
                0,sy,0;... 4
                0,0,sz;... 5
                sx,0,sz;... 6
                sx,sy,sz;... 7
                0,sy,sz]; % 8
            faces=[1,2;...
                2,3;...
                1,4;...
                3,4;...
                1,5;...
                5,6;...
                5,8;...
                6,7;...
                7,8;...
                2,6;...
                3,7;...
                4,8];
        end
        function [nodes, faces]= bccFcc(sx, sy, sz)
            nodes=[0,0,0;...
                (sx),0,0;...
                sx,sy,0;...
                0,sy,0;0,0,sz;...
                sx,0,sz;...
                sx,sy,sz;...
                0,sy,sz;...
                sx/2,sy/2,sz/2;...
                sx/2,sy/2,0;...
                sx/2,sy/2,sz;...
                sx/2,0,sz/2;...
                sx/2,sy,sz/2;...
                0,sy/2,osz/2;...
                sx,sy/2,sz/2];
        end
        function [nodes, faces]= fbc(sx, sy, sz)
            nodes=[0,0,0;...
                sx,0,0;...
                0,sy,0;...
                sx,sy,0;...
                0,0,sz;...
                sx,0,sz;...
                0,sy,sz;...
                sx,sy,sz];
            faces=[2,8;4,6;1,7;3,5;1,6;2,5;4,7;3,8;... faces
                1,8;2,7;4,5;3,6]; %body
        end
        function [nodes, faces]= fbcxyz(sx, sy, sz)
            nodes=[0,0,0;...
                sx,0,0;...
                0,sy,0;...
                sx,sy,0 ;...
                0,0,sz;...
                sx,0 ,sz;...
                0,sy,sz;...
                sx,sy,sz];
            faces=[1,2;1,3;4,3;4,2;2,6;6,8;8,4;1,5;5,7;7,3;5,6;7,8;... edges
                1,4;2,3;2,8;4,6;1,7;3,5;1,6;2,5;4,7;3,8;7,6;5,8;... faces
                1,8;2,7;4,5;3,6]; %body
        end
        function [nodes, faces]= fbcz(sx, sy, sz)
            nodes=[0,0,0;...
                sx,0,0;...
                0,sy,0;...
                sx,sy,0;...
                0,0,sz;...
                sx,0,sz;...
                0,sy,sz;...
                sx,sy,sz];
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
    methods (Static) % strut splitting methods
        function boundingBox = getBound(vertices,faces,desired)
            % get bounding box of ervery strut
            %[minX minY minZ]
            %[maxX maxY maxZ]
            boundingBox = cell(length(desired),1);
            numRequiredBB = length(desired);
            for inc = 1:numRequiredBB
                index = desired(inc);
                maxers = max(vertices(faces(index,1),:),vertices(faces(index,2),:));
                miners = min(vertices(faces(index,1),:),vertices(faces(index,2),:));
                boundingBox{inc} = [miners;maxers];
            end
        end
        function isIntersect = findBbIntersect(currentBB,boundingBox)
            % determine which boundingBox intersect the currentBB
            numBB = length(boundingBox);
            isIntersect = false(numBB,1);
            for inc = 1:numBB
                checkBB = boundingBox{inc};
                %check each dimension
                inX = PLG.checkIn(currentBB(:,1),checkBB(:,1));
                inY = PLG.checkIn(currentBB(:,2),checkBB(:,2));
                inZ = PLG.checkIn(currentBB(:,3),checkBB(:,3));
                
                isIntersect(inc) = inX & inY & inZ;
            end
        end
        function potentialStruts = removeNormalConStruts(p0,p1,potentialStruts,faces,verts,tol)
            % remove struts that are already joined properly
            numPotential = length(potentialStruts);
            test = false(numPotential,1);
            for inc = 1:numPotential
                index = potentialStruts(inc);
                q0 = verts(faces(index,1),:);
                q1 = verts(faces(index,2),:);
                diff1 = all(abs(p0-q0)<tol);
                diff2 = all(abs(p1-q0)<tol);
                diff3 = all(abs(p0-q1)<tol);
                diff4 = all(abs(p1-q1)<tol);
                test(inc) = any([diff1,diff2,diff3,diff4]);
            end
            potentialStruts(test) = [];
        end
        function potentialStruts = removeParralelStruts(p0,p1,potentialStruts,faces,verts)
            v = (p0-p1)/norm(p0-p1);
            numPotential = length(potentialStruts);
            test = false(numPotential,1);
            for inc = 1:numPotential
                index = potentialStruts(inc);
                q0 = verts(faces(index,1),:);
                q1 = verts(faces(index,2),:);
                u = (q0-q1)/norm(q0-q1);
                test(inc) = all(abs(v-u)<1e-6);
            end
            potentialStruts(test) = [];
        end
        function potentialStruts =  findSplitStrutBB(splitStruts,potentialStruts,currentBB,currentP1,currentP2,facesOut,vertsOut,tol)
            % return split strut and orriginal potenial intersecting struts
            potentialStruts = potentialStruts(:);
            numPotentialStruts = length(potentialStruts);
            potentialStruts = [potentialStruts,ones(numPotentialStruts,1)]; % second column will be 1 for orriginal array 2 for output array
            % determine if any of the above potential struts have already being split
            test = ~cellfun(@isempty,splitStruts(potentialStruts(:,1)));
            strutsToReplace = potentialStruts(test,1);
            
            if isempty(strutsToReplace)
                % no changes needed
            else
                % replace whole struts with split ones
                % remove whole struts that are to be replaced
                potentialStruts(test,:) = [];
                % get split struts BB
                desired = [];
                for inc = strutsToReplace'
                    desired = [desired,splitStruts{inc}(1):splitStruts{inc}(2)];
                end
                boundingBox = PLG.getBound(vertsOut,facesOut,desired);
                isIntersect = PLG.findBbIntersect(currentBB,boundingBox);
                potentialNewStruts = desired(isIntersect);
                potentialNewStruts = PLG.removeNormalConStruts(currentP1,currentP2,potentialNewStruts,facesOut,vertsOut,tol);
                % dont need to remove parallel as these are already removed in whole strut check
                if isempty(potentialNewStruts)
                    % nothing
                else
                    potentialStruts = [potentialStruts;...
                        [potentialNewStruts',ones(length(potentialNewStruts),1)*2]];
                end
                
            end
        end
        function [vertsTmpOut,facesTmpOut]= findIntersect(p0,p1,potentialStruts,faces,verts,facesOut,vertsOut,tol)
            % for struts that have a coincident bounding box determine if they
            % intersect
            numPotentialStruts = size(potentialStruts,1);
            coordinates = zeros(numPotentialStruts,3);
            order = zeros(numPotentialStruts,1);
            for inc = 1:numPotentialStruts
                strutType = potentialStruts(inc,2);
                index = potentialStruts(inc,1);
                if strutType==1
                    q0 = verts(faces(index,1),:);
                    q1 = verts(faces(index,2),:);
                else
                    q0 = vertsOut(facesOut(index,1),:);
                    q1 = vertsOut(facesOut(index,2),:);
                end
                [dist,sc,p1Out,p2Out] = PLG.distBetween2Segment(p0, p1, q0, q1);
                if dist<tol
                    coordinates(inc,:) = mean([p1Out;p2Out],1);
                    order(inc) = sc; % the larger the value the further the point is from p1
                else
                    coordinates(inc,:) = NaN;
                    order(inc) = NaN; % the larger the value the further the point is from p1
                end
            end
            % get unique coordinates with NaN removed
            test = isnan(order);
            order(test)=[];
            coordinates(test,:) = [];
            [order,index] = uniquetol(order,tol,'DataScale',1);
            coordinates = coordinates(index,:);
            % create a face and verts array
            numNewFaces = length(order)+1;
            vertsTmpOut = [p1;coordinates;p0];
            facesTmpOut = [1:numNewFaces;2:(numNewFaces+1)]';
        end
        function isIn = checkIn(pointsMain,pointsCheck)
            % determine if pointsCheck intersect pointsMain
            % pointsMain = [minPoint,maxPoint]
            if pointsMain(1) < pointsCheck(1)
                lowRange = pointsMain;
                highRange = pointsCheck;
            else
                lowRange = pointsCheck;
                highRange = pointsMain;
            end
            
            isIn = lowRange(2) >= highRange(1);
            
        end
        function [distance fromP1 coord1 coord2] = distBetween2Segment(p0, p1, q0, q1)
            % code from: https://au.mathworks.com/matlabcentral/fileexchange/32487-shortest-distance-between-two-line-segments?focused=3821416&tab=function
            u = p0 - p1;
            v = q0 - q1;
            w = p1 - q1;
            
            a = dot(u,u);
            b = dot(u,v);
            c = dot(v,v);
            d = dot(u,w);
            e = dot(v,w);
            D = a*c - b*b;
            sD = D;
            tD = D;
            
            SMALL_NUM = 0.00000001;
            
            % compute the line parameters of the two closest points
            if (D < SMALL_NUM)  % the lines are almost parallel
                sN = 0.0;       % force using point P0 on segment S1
                sD = 1.0;       % to prevent possible division by 0.0 later
                tN = e;
                tD = c;
            else                % get the closest points on the infinite lines
                sN = (b*e - c*d);
                tN = (a*e - b*d);
                if (sN < 0.0)   % sc < 0 => the s=0 edge is visible
                    sN = 0.0;
                    tN = e;
                    tD = c;
                elseif (sN > sD)% sc > 1 => the s=1 edge is visible
                    sN = sD;
                    tN = e + b;
                    tD = c;
                end
            end
            
            if (tN < 0.0)            % tc < 0 => the t=0 edge is visible
                tN = 0.0;
                % recompute sc for this edge
                if (-d < 0.0)
                    sN = 0.0;
                elseif (-d > a)
                    sN = sD;
                else
                    sN = -d;
                    sD = a;
                end
            elseif (tN > tD)       % tc > 1 => the t=1 edge is visible
                tN = tD;
                % recompute sc for this edge
                if ((-d + b) < 0.0)
                    sN = 0;
                elseif ((-d + b) > a)
                    sN = sD;
                else
                    sN = (-d + b);
                    sD = a;
                end
            end
            
            % finally do the division to get sc and tc
            if(abs(sN) < SMALL_NUM)
                sc = 0.0;
            else
                sc = sN / sD;
            end
            
            if(abs(tN) < SMALL_NUM)
                tc = 0.0;
            else
                tc = tN / tD;
            end
            
            % get the difference of the two closest points
            dP = w + (sc * u) - (tc * v);  % = S1(sc) - S2(tc)
            
            %% outputs
            distance = norm(dP);
            fromP1 = sc;
            coord1 = p1+sc*u;   % Closest point on object 1
            coord2 = q1+tc*v;   % Closest point on object 2
            
        end
    end
end
