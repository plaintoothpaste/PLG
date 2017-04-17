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
    properties (SetAccess=private)
        latticeType;
        resolution;
        strutDiamter;
        % sphereAddition;
        sphereDiameter;
        usx;
        usy;
        usz;
        repsx;
        repsy;
        repsz;
        orginx;
        orginy;
        orginz;
        
        latticeStructure;
    end
    properties (Access = protected)
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
                    
                    obj.strutureType = 1;
                otherwise
                    error('Incorrect number of inputs');
            end
        end
        function obj = latticeGenerate(obj)
            % generates a latticeStructure ready for simulation or saving as an stl etc
            origin = [obj.orginx,obj.orginy,obj.orginz];
            switch obj.latticeType
                case 'bcc' % BCC Cell
                    [vertices, faces] = PLG.bcc(obj.usx,obj.usy,obj.usz,origin);
                case 'bcc2' % BCC Cell with no cantilevers
                    [vertices, faces] = PLG.bcc2(obj.usx,obj.usy,obj.usz,origin);
                case 'bcz' % BCC Cell
                    [vertices, faces] = PLG.bcz(obj.usx,obj.usy,obj.usz,origin);
                case 'fcc'
                    [vertices, faces] = PLG.fcc(obj.usx,obj.usy,obj.usz,origin);
                case 'fccNoXY'
                    [vertices, faces] = PLG.fccNoXY(obj.usx,obj.usy,obj.usz,origin);
                case 'fbcxyz'
                    [vertices, faces] = PLG.fbcxyz(obj.usx,obj.usy,obj.usz,origin);
                case 'fcz'
                    [vertices, faces] = PLG.fcz(obj.usx,obj.usy,obj.usz,origin);
                case 'fbcz'
                    [vertices, faces] = PLG.fbcz(obj.usx,obj.usy,obj.usz,origin);
                case 'bcc_fcc'
                    [vertices, faces] = PLG.bccFcc(obj.usx,obj.usy,obj.usz,origin);
                case 'fbc'
                    [vertices, faces] = PLG.fbc(obj.usx,obj.usy,obj.usz,origin);
                case 'box'
                    [vertices, faces] = PLG.box(obj.usx,obj.usy,obj.usz,origin);
            end
            % replicate the unit cells generated above
            [obj.latticeStructure.vertices, obj.latticeStructure.faces] = cellReplication(obj,vertices,faces);
            
            
        end
        function [vertices, faces] = cellReplication(obj,vertices,faces)
            count = 1;
            no_nodes=size(vertices,1);
            node_connect=size(faces,1);
            no_links=max(faces);
            no_links=max(no_links);
            vert_blank=zeros(obj.repsx*obj.repsy*obj.repsz*(no_nodes),3); %prealocate arrays
            face_blank=zeros(obj.repsx*obj.repsy*obj.repsz*node_connect,2); %prealocate arrays
            vert_add=zeros(no_nodes,3);
            
            %completion notification
            for i=0:obj.repsx-1
                %completion notification
                vert_add(:,1)=vertices(:,1)+i*obj.usx;
                for j=0:obj.repsy-1
                    vert_add(:,2)=vertices(:,2)+j*obj.usy;
                    for k=0:obj.repsz-1
                        vert_add(:,3)=vertices(:,3)+k*obj.usz;
                        vert_blank((count-1)*no_nodes+1:((count)*no_nodes),:)=vert_add;
                        face_blank((count-1)*node_connect+1:((count)*node_connect),:)=faces+(count-1)*(no_links);
                        count=count+1;
                    end
                end
            end
            % remove duplicate verts and clean up
            [vertices,~,indexn]=unique(vert_blank,'rows');
            faces = indexn(face_blank);
            
            % remove duplicate faces
            for inc = 1:length(faces)
                ind1 = faces(inc,1);
                ind2 = faces(inc,2);
                if ind1<=ind2
                    % faces(inc,1) = ind1;
                    % faces(inc,2) = ind2;
                    % no change
                else
                    faces(inc,1) = ind2;
                    faces(inc,2) = ind1;
                end
            end
            faces = unique(faces,'rows');
        end
        function obj = loadIn(obj,file)
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
                    obj.latticeStructure.vertices = data(3:numNodes+2,1:3);
                    
                    obj.latticeStructure.faces    = data(numNodes+3:numNodes+numLinks+2,1:2);
                    obj.latticeStructure.diameters = data(numNodes+3:numNodes+numLinks+2,3);
                    if size(data,2)==3
                        % no sphere diameter supplied
                        obj.latticeStructure.spheres = zeros(numNodes,1);
                    else
                        % sphere diameter supplied
                        obj.latticeStructure.spheres = data(3:numNodes+2,4);
                    end
                case obj.validExtensions{3}
                    % custom
                    data = csvread(file);
                    numNodes=data(1,1);
                    numLinks=data(2,1);
                    % data
                    obj.latticeStructure.vertices = data(3:numNodes+2,1:3);
                    obj.latticeStructure.faces    = data(numNodes+3:numNodes+numLinks+2,1:2);
                    obj.latticeStructure.diameters = data(numNodes+3:numNodes+numLinks+2,3);
                    if size(data,2)==3
                        % no sphere diameter supplied
                        obj.latticeStructure.spheres = zeros(numNodes,1);
                    else
                        % sphere diameter supplied
                        obj.latticeStructure.spheres = data(3:numNodes+2,4);
                    end
            end
        end
        function obj = translate(obj,x,y,z)
            obj.latticeStructure.vertices(:,1) = obj.latticeStructure.vertices(:,1)+x;
            obj.latticeStructure.vertices(:,2) = obj.latticeStructure.vertices(:,2)+y;
            obj.latticeStructure.vertices(:,3) = obj.latticeStructure.vertices(:,3)+z;
        end
        function obj = cleanLattice(obj)
            % cleans the lattice structure including the removal of
            % duplicate vertices and struts
            % remove duplicate faces
            for inc = 1:length(obj.latticeStructure.faces)
                ind1 = obj.latticeStructure.faces(inc,1);
                ind2 = obj.latticeStructure.faces(inc,2);
                if ind1<=ind2
                    % faces(inc,1) = ind1;
                    % faces(inc,2) = ind2;
                    % no change
                else
                    obj.latticeStructure.faces(inc,1) = ind2;
                    obj.latticeStructure.faces(inc,2) = ind1;
                end
            end
            [obj.latticeStructure.faces,i,j] = unique(obj.latticeStructure.faces,'rows');
            obj.latticeStructure.diameters = obj.latticeStructure.diameters(i);
            
            [obj.latticeStructure.vertices,i,indexn]=unique(obj.latticeStructure.vertices,'rows');
            obj.latticeStructure.spheres = obj.latticeStructure.spheres(i);
            obj.latticeStructure.faces = indexn(obj.latticeStructure.faces);
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
            numPoints = length(obj.latticeStructure.vertices);
            newPoints = zeros(size(obj.latticeStructure.vertices));
            for inc = 1:numPoints
                newPoints(inc,:) = obj.latticeStructure.vertices(inc,:)*rx';
            end
            for inc = 1:numPoints
                
            end
            for inc = 1:numPoints
                
            end
            obj.latticeStructure.vertices = newPoints;
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
            
            PLG.plotLine(obj.latticeStructure.vertices,obj.latticeStructure.faces);
            PLG.plotPoints(obj.latticeStructure.vertices);
            
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
                    % TODO
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
            fullName = [pathName,fileName];
            numFacets=obj.resolution*4;%number of facets created for one strut
            numLinks=size(obj.latticeStructure.faces,1);
            numVertices=size(obj.latticeStructure.vertices,1);
            strutFacets=numFacets*numLinks; %total number of facets for all the struts
            
            fid=fopen(fullName,'w');
            fprintf(fid, '%-80s', 'fast stl generator'); %binary write information
            
            if ~isempty(obj.sphereDiameter)
                fwrite(fid,uint32(strutFacets+8^2*(8-1)*numVertices),'uint32'); %stl binary header file contains the total number of facets in the stl file
                fid = PLG.ballCreate(obj.latticeStructure,obj.sphereDiameter/2,numVertices,fid);
            else
                fwrite(fid,uint32(strutFacets),'uint32'); %stl binary header file contains the total number of facets in the stl file
            end
            fid = PLG.faceCreate(obj.latticeStructure,obj.resolution,obj.strutDiamter/2,numLinks,fid);
            fclose(fid);
        end
        function saveExcel(obj,fileName,pathName)
            % saves data to an excel compatible format
            fullName = [pathName,fileName];
            
            numNodes=size(obj.latticeStructure.vertices,1);
            numLinks=size(obj.latticeStructure.faces,1);
            xlswrite(fullName,numNodes,'Sheet1','A1');
            xlswrite(fullName,numLinks,'Sheet1','A2');
            
            locationVertices=strcat('A3:C',num2str(numNodes+2));
            xlswrite(fullName,obj.latticeStructure.vertices,'Sheet1',locationVertices);
            locationFaces=strcat('A',num2str(numNodes+3),':B',num2str(numNodes+2+numLinks));
            xlswrite(fullName,obj.latticeStructure.faces,'Sheet1',locationFaces);
            
            % optional write depending on load type
            if obj.strutureType ==0
                % loaded custom
                locationSpheres = strcat('D3:D',num2str(numNodes+2));
                xlswrite(fullName,obj.latticeStructure.spheres,'Sheet1',locationSpheres);
                locationDiameters = strcat('C',num2str(numNodes+3),':C',num2str(numNodes+2+numLinks));
                xlswrite(fullName,obj.latticeStructure.diameters,'Sheet1',locationDiameters);
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
            
            numNodes=size(obj.latticeStructure.vertices,1);
            numLinks=size(obj.latticeStructure.faces,1);
            
            dlmwrite(fullName,numNodes);
            dlmwrite(fullName,numLinks,'-append');
            data = [obj.latticeStructure.vertices,obj.latticeStructure.spheres];
            dlmwrite(fullName,data,'-append');
            
            data = [obj.latticeStructure.faces,obj.latticeStructure.diameters];
            dlmwrite(fullName,data,'-append');
        end
    end
    methods (Static) % non cell type methods
        function  fid = ballCreate(latticeStructure,radius,numVertices,fid)
            [x,y,z]=sphere(8); %create sphere with higher accuracy
            x=x(:);y=y(:);z=z(:); %reshape points
            ball.faces= convhull([x y z]); %create triangle links
            sizer = size(ball.faces,1);
            ball.vertices=[x,y,z]; %store the points
            offset=ball.vertices*radius;
            for i=1:numVertices
                target=[offset(:,1)+latticeStructure.vertices(i,1),offset(:,2)+latticeStructure.vertices(i,2),offset(:,3)+latticeStructure.vertices(i,3)];
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
        function fid = faceCreate(latticeStructure,resolution,radius,numLinks,fid)
            numVertices=((resolution+1)*2); %number of points created for 1 strut
            
            %% create facets linking information
            % for face 1 create an array holding all the facet links for a strut's end
            end1=ones(resolution,3);
            end1(:,2)=[2:resolution+1]';
            end1(:,3)=[3:resolution+1,2]';
            % for face 2
            end2=end1+resolution+1;
            %sides in both directions
            end1end2=[end1(:,2),end2(:,2:3)];
            end2end1=[end1(:,2:3),end2(:,3)];
            
            for i=1:numLinks
                point1=latticeStructure.vertices(latticeStructure.faces(i,1),:);
                point2=latticeStructure.vertices(latticeStructure.faces(i,2),:);
                separation=(point2-point1); %point1's vector normalised
                u1=separation/norm(separation);
                %create offset vector
                if u1(3)~=0
                    Vs=cross([1,0,0],u1);
                    v1=Vs/norm(Vs);
                elseif u1(2)~=0
                    Vs=cross([0,0,1],u1);
                    v1=Vs/norm(Vs);
                elseif u1(1)~=0
                    Vs=cross([0,1,0],u1);
                    v1=Vs/norm(Vs);
                else
                    error('the two nodes collide')
                end
                offset1=radius*cross(v1,u1)';
                for j=1:(resolution+1)*2
                    if j==1
                        vertices(j,:)=point1;
                        count=1;
                    elseif j>1 && j<=resolution+1
                        Qrot1 = PLG.qGetRotQuaternion( ((count*2*pi)/resolution), u1 );
                        vertices(j,:)=PLG.qRotatePoint( offset1 , Qrot1 )'+point1;
                        
                        %vertices(j,:) = rad*cos(alpha*count)*v1+rad*sin(alpha*count)*cross(v1,u1)+point1;
                        count=count+1;
                    elseif j>resolution+1
                        vertices(j,:)=vertices(j-(resolution+1),:)+separation;
                    end
                end
                
                %% join faces
                for j=1:resolution
                    %get values first end
                    facet_a=vertices(end1(j,3),:);
                    facet_b=vertices(end1(j,2),:);
                    facet_c=vertices(end1(j,1),:);
                    normal=cross(facet_b-facet_a,facet_c-facet_a);
                    %write values
                    fwrite(fid,normal,'float32');           % write normal vector floating point numbers
                    fwrite(fid,facet_a,'float32');   % first vertex (x,y,z)
                    fwrite(fid,facet_b,'float32');   % second vertex
                    fwrite(fid,facet_c,'float32');   % Third vertex
                    fwrite(fid,0,'uint16','l');
                    
                    %get values second end
                    facet_a=vertices(end2(j,1),:);
                    facet_b=vertices(end2(j,2),:);
                    facet_c=vertices(end2(j,3),:);
                    normal=cross(facet_b-facet_a,facet_c-facet_a);
                    %write values
                    fwrite(fid,normal,'float32');           % write normal vector floating point numbers
                    fwrite(fid,facet_a,'float32');   % first vertex (x,y,z)
                    fwrite(fid,facet_b,'float32');   % second vertex
                    fwrite(fid,facet_c,'float32');   % Third vertex
                    fwrite(fid,0,'uint16','l');
                    
                    %get values cylinder 1 direction
                    facet_a=vertices(end1end2(j,3),:);
                    facet_b=vertices(end1end2(j,2),:);
                    facet_c=vertices(end1end2(j,1),:);
                    normal=cross(facet_b-facet_a,facet_c-facet_a);
                    %write values
                    fwrite(fid,normal,'float32');           % write normal vector floating point numbers
                    fwrite(fid,facet_a,'float32');   % first vertex (x,y,z)
                    fwrite(fid,facet_b,'float32');   % second vertex
                    fwrite(fid,facet_c,'float32');   % Third vertex
                    fwrite(fid,0,'uint16','l');
                    
                    %get values cylinder in other direction
                    facet_a=vertices(end2end1(j,1),:);
                    facet_b=vertices(end2end1(j,2),:);
                    facet_c=vertices(end2end1(j,3),:);
                    normal=cross(facet_b-facet_a,facet_c-facet_a);
                    %write values
                    fwrite(fid,normal,'float32');           % write normal vector floating point numbers
                    fwrite(fid,facet_a,'float32');   % first vertex (x,y,z)
                    fwrite(fid,facet_b,'float32');   % second vertex
                    fwrite(fid,facet_c,'float32');   % Third vertex
                    fwrite(fid,0,'uint16','l');
                end
            end
        end
        function plotLine(vert,strut)
            % plot a heap of lines fast
            numStrut = size(strut,1);
            for inc = 1:numStrut
                link = strut(inc,:);
                x = [vert(link(1),1),vert(link(2),1)];
                y = [vert(link(1),2),vert(link(2),2)];
                z = [vert(link(1),3),vert(link(2),3)];
                l = line(x,y,z);
                l.Color	=[0.3,0.3,0.3,0.5];
            end
        end
        function plotPoints(vert)
            % plot a heap of verts
            if isempty(vert) % edge
                return
            end
            numVert = size(vert,1);
            for inc = 1:numVert
                points = vert(inc,:);
                s = scatter3(points(1),points(2),points(3));
                s.MarkerFaceColor = [0.9,0.5,0];
                s.MarkerEdgeColor = 'none';
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
            nodes=[origin;...                                           1  0
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
