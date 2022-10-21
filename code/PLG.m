classdef PLG
    %PLG a program designed for the generation of lattice structures.
    %    The PLG program designed for the generation of regular repeating 
    %    stuctures, normally a lattice. See readme for details.
    
    properties (SetAccess=protected)
        unitName;
        
        resolution;
        strutDiameter;
        
        sphereAddition;
        sphereResolution;
        sphereDiameter;
        baseFlat; % logical value that only matters if sphere addition is also true, deforms sphere to have a flat base at minZ
        
        vertices;
        struts;
        
        % 3 value vectors for the x y z params respectivally
        unitSize
        replications
        origin
    end
    properties (SetAccess=private)
        runLocation; % the location from which the PLG is running.
        
        loadExtensions = {'xml','custom unit cell file defined as a xml';...
                          'stl', 'standard AM triangular surface file';...
                          'lattice','standard beam output of PLG';...
                          'csv','manualy defined beam model';...
                          'xlsx','manualy defined beam model'};
        saveExtensions = {'*.stl', 'binary facet representation (compatibility)';...
                          '*.3mf', 'tesselated facet file (recommended)';...
                          '*.lattice', 'lattice beam output method (to use in PLG in the future)';...
                          '*.inp', 'Abaqus input file'};
        
        tolerance; % defined as 1/100 of the shortest length present
    end
    methods
        function obj = PLG(varargin)
            % PLG generates a lattice or loads an existing one.
            %   With no inputs a new lattice is made
            %   With one input of a file name the '.lattice' file is loaded
            %   Initialises the PLG object
            switch numel(varargin)
                case 0
                    % generate a new lattice
                case 1
                    % import a custom lattice file containing beam and node
                    % definitions see load function for more information
                    obj = load(obj,varargin{1});
                    
                    obj.unitName = varargin{1};
                otherwise
                    error('Only a single file path can be specified.');
            end
            obj.sphereAddition = false;
            obj.baseFlat = false;
            obj.origin=[0,0,0];
            
            pathPLG = which('PLG');
            obj.runLocation = fileparts(pathPLG);
        end
        function obj = set(obj,varargin)
            % set set properties of the PLG object
            %    contains some basic validation.
            %    Returns a modified PLG object.
            allowable={'resolution','strutDiameter','sphereAddition',...
                'sphereResolution','sphereDiameter','unitSize',...
                'replications','origin','baseFlat'};
            classType = {'double','double','logical',...
                'double','double','double',...
                'double','double','logical'};
            attribute = {{'scalar','nonempty'},{'scalar','nonempty'},{'scalar','nonempty'},...
                {'scalar','nonempty'},{'scalar','nonempty'},{'size',[1,3],'nonempty'},...
                {'size',[1,3],'nonempty'},{'size',[1,3],'nonempty'},{'scalar','nonempty'}};
            % setup parser
            p = inputParser();
            for inc = 1:length(allowable)
                name = allowable{inc};
                class = classType{inc};
                attr = attribute{inc};
                f = @(x) validateattributes(x,{class},attr);
                addParameter(p,name,[],f)
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
        function obj = defineUnit(obj,unitNames)
            % defineUnit creates the unit cell that is then manipulated by PLG
            %    unitNames - A list of strings with each desired unit cell.
            %    Returns a modified PLG object.
            %
            %    WARNING: Running this function will delete any replication
            %    translation rotation etc that were were previously applied.
            addpath([obj.runLocation,filesep,'unitCell']);
            unitObj = unitCell(unitNames,obj);
            rmpath([obj.runLocation,filesep,'unitCell']);
            
            obj.vertices = unitObj.vertices;
            obj.struts = unitObj.connections;
            obj.strutDiameter = unitObj.strutDiam;
            obj.sphereDiameter = unitObj.sphereDiam;
            obj.unitName = unitObj.unitName;
        end
        function obj = cellReplication(obj)
            % cellReplication Applies cell replications.
            %    The details of cell replication must be set before calling
            %    this function.
            %    Returns a modified PLG object.
            
            xPlacement = [0:obj.unitSize(1):obj.unitSize(1)*(obj.replications(1)-1)]+obj.origin(1);
            yPlacement = [0:obj.unitSize(2):obj.unitSize(2)*(obj.replications(2)-1)]+obj.origin(2);
            zPlacement = [0:obj.unitSize(3):obj.unitSize(3)*(obj.replications(3)-1)]+obj.origin(3);
            [XX,YY,ZZ] = ndgrid(xPlacement,yPlacement,zPlacement);
            
            vertOut = arrayfun(@(x,y,z) ...
                [obj.vertices(:,1) + x,...
                obj.vertices(:,2) + y,...
                obj.vertices(:,3) + z]...
                ,XX,YY,ZZ,'UniformOutput',0);
            obj.vertices = cell2mat(vertOut(:));
            
            sphereDiamOut = arrayfun(@(x) ...
                obj.sphereDiameter,XX,...
                'UniformOutput',0);
            obj.sphereDiameter = cell2mat(sphereDiamOut(:));
            
            
            strutCounter = 0:(numel(XX)-1);
            numStruts = max(max(obj.struts));
            strutCounter = strutCounter*numStruts;
            strutOut = arrayfun(@(x) obj.struts + x,strutCounter,'UniformOutput',0);
            obj.struts = cell2mat(strutOut(:));
            strutDiamOut = arrayfun(@(x) obj.strutDiameter, XX,'UniformOutput',0);
            obj.strutDiameter = cell2mat(strutDiamOut(:));
        end
        function obj = cleanLattice(obj,tol)
            % cleanLattice removes coincident vertices and struts.
            %    Due to transformations, replications etc there may be
            %    items that coincide with each other. This function
            %    removes those items. It is recommended to always run this
            %    function before saving out a file.
            %    Returns a modified PLG object.
            %
            %    WARNING: may have a strange effect on diameters.
            
            verts1 = obj.vertices(obj.struts(:,1),:);
            verts2 = obj.vertices(obj.struts(:,2),:);
            lengthVerts = sum(sqrt((verts1-verts2).^2),2);
            
            if ~exist('tol','var')
                obj.tolerance = min(lengthVerts)/20;
            else
                obj.tolerance = tol;
            end
            
            
            % duplicate vertices
            [obj.vertices,i,indexn]=uniquetol(obj.vertices,obj.tolerance,'ByRows',1,'DataScale',1);
            obj.sphereDiameter = obj.sphereDiameter(i);
            obj.struts = indexn(obj.struts);
            
            % duplicate struts
            tmp = obj.struts;
            test = tmp(:,1)>tmp(:,2);
            obj.struts(test,1) = tmp(test,2);
            obj.struts(test,2) = tmp(test,1);
            [obj.struts,i] = unique(obj.struts,'rows');
            obj.strutDiameter = obj.strutDiameter(i);
            
            % duplicate struts zero length
            test = obj.struts(:,1)==obj.struts(:,2);
            obj.struts(test,:)=[];
            obj.strutDiameter(test)=[];
            
            % remove unused verts
            numVerts = size(obj.vertices,1);
            allInds = (1:numVerts)';
            usedInds = unique(obj.struts(:));
            unusedInds = setdiff(allInds,usedInds);
            if ~isempty(unusedInds)
                newInds = (1:length(usedInds))';
                allInds(usedInds) = newInds;
                allInds(unusedInds)=0;
                obj.struts = allInds(obj.struts);
                obj.vertices(unusedInds,:) = [];
                obj.sphereDiameter(unusedInds) = [];
            end
        end
        function [f,a] = plot(obj,colours)
            % plot makes a basic plot of the lattices 
            %    colours - a single colour for all struts or a nx3 array of
            %    colours where n is the length of struts.
            %    Returns f-figure handle.
            %    Returns a-axes handle.
            f = figure;
            f.Units	= 'Normalized';
            f.Position = [0.1,0.1,0.8,0.8];
            f.Name = 'PLG plotter';
            a = axes;
            a.View = [45,45];
            axis vis3d
            a.NextPlot='add';
            numStruts = length(obj.strutDiameter);
            if ~exist('colours','var')
                colours = repmat([0.1,0.1,1.0,0.5],numStruts,1);
            end
            for inc = 1:numStruts
                points = [obj.vertices(obj.struts(inc,1),:);obj.vertices(obj.struts(inc,2),:)];
                p = plot3(points(:,1),points(:,2),points(:,3),'Color',colours(inc,:));
                p.MarkerFaceColor = [0.9,0.5,0]; p.MarkerEdgeColor = 'none';
                p.Marker = 'o';
                p.LineWidth = 3;
            end
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
        function save(obj)
            % this overloads the save function and allows saving out in various
            % formats however this only saves the latticeStructure structure
            [fileName,pathName,filterIndex] = uiputfile(obj.saveExtensions);
            file = [pathName,fileName];
            switch filterIndex
                case 1
                    saveStl(obj,file);
                case 2
                    save3mf(obj,file)
                case 3
                    saveLattice(obj,file)
                case 4
                    saveAbaqus(obj,file)
                otherwise
                    error('No file saved')
            end
        end
    end
    methods % lattice manipulations
        function obj = scale(obj,sx,sy,sz)
            % scale Scales the struts must provide a scale in x,y,z
            %    Does not scale diameters.
            %    Returns a modified PLG object.
            obj.vertices(:,1) = obj.vertices(:,1)*sx;
            obj.vertices(:,2) = obj.vertices(:,2)*sy;
            obj.vertices(:,3) = obj.vertices(:,3)*sz;
        end
        function obj = translate(obj,x,y,z)
            % translate Translate the struts must provide a x,y,z
            %    Returns a modified PLG object.
            obj.vertices(:,1) = obj.vertices(:,1)+x;
            obj.vertices(:,2) = obj.vertices(:,2)+y;
            obj.vertices(:,3) = obj.vertices(:,3)+z;
        end
        function obj = rotate(obj,wx,wy,wz)
            % rotate Rotate the struts must provide a rotation about x,y,z
            %    Rotations are in degrees about the main axes
            %    Returns a modified PLG object.
            
            % to radians
            wx = wx*pi/180;
            wy = wy*pi/180;
            wz = wz*pi/180;
            
            affineMatrix = [1*cos(wy)*cos(wz), sin(wz),           -sin(wy),                 0;...
                            -sin(wz),          cos(wx)*1*cos(wz), sin(wx),           0;...
                            sin(wy),           -sin(wx),          cos(wx)*cos(wy)*1, 0;...
                            0,                 0,                 0,                 1]; % needs to be transposed to work with row rather then column data that the verts are stored in 
            newVerts = arrayfun(@(x,y,z) [x,y,z,1]*affineMatrix,obj.vertices(:,1),obj.vertices(:,2),obj.vertices(:,3),'UniformOutput',0);
            newVerts = cell2mat(newVerts);
            obj.vertices = newVerts(:,1:3);
        end
        function obj = plus(obj,obj1)
            % plus Combines two objects together. Overloads the "+"
            %    Returns a modified PLG object.
            %
            %    WARNING: may have a strange effect on diameters.
            newStruts = obj1.struts+size(obj.vertices,1); % in case there are some verts not being used
            obj.struts = [obj.struts;newStruts];
            obj.strutDiameter = [obj.strutDiameter;obj1.strutDiameter];
            obj.sphereDiameter = [obj.sphereDiameter;obj1.sphereDiameter];
            obj.vertices = [obj.vertices;obj1.vertices];
            
            % remove any matching struts
            obj = cleanLattice(obj);
        end
    end
    methods % stats and advanced
        function obj = calcDx(obj)
            % calcDx Calculate the length of all struts.
            p1 = obj.vertices(obj.struts(:,1),:);
            p2 = obj.vertices(obj.struts(:,2),:);
            diffx = (p1(:,1)-p2(:,1));
            diffy = (p1(:,2)-p2(:,2));
            diffz = (p1(:,3)-p2(:,3));
            % TODO
            obj.dx = sqrt(diffx.^2+diffy.^2+diffz.^2);
        end
        function getProperties(obj,fileName)
            % getProperties Perform a statistical analysis.
            % convert plg inputs into the function inputs
            % TODO move into a sub class
            numNodes=size(obj.vertices,1);
            numStruts=size(obj.struts,1); % Read Number of Struts
            Struts=[obj.struts,obj.strutDiameter];
            Nodes=[obj.vertices(:,1) obj.vertices(:,3) obj.vertices(:,2) obj.sphereDiameter];
            
            %% CHECK LATTICE DIMENSIONALITY
            % Checks if nodes lie in a single plane
            % 1=3D and 2=2D
            if length(unique(Nodes(:,1)))==1 || length(unique(Nodes(:,2)))==1 || length(unique(Nodes(:,3)))==1
                Dimen=2;
                M=numStruts-(2*numNodes)+3;
            else
                Dimen=1;
                M=numStruts-(3*numNodes)+ 6;
            end
            % Calculate Maxwell number and categorise
            if M>0
                LatBehav='Over-stiff';
            elseif M==0
                LatBehav='Just-stiff';
            elseif M<0
                LatBehav='Under-stiff';
            end
            % Write for Output Excel Sheet format
            LatBehavOut=sprintf('%d (%s)',M,LatBehav);
            
            %% Analyse Struts
            % Create Counters for Strut Manufacture Quality
            RobustZoneCounter=0;
            CompromisedZoneCounter=0;
            FailedZoneCounter=0;
            summaryData=cell(0);
            for idxstrut=1:numStruts
                %%CURRENT STRUT DETAILS
                curStrut=Struts(idxstrut,:);
                StrutDia=curStrut(3); % Current strut diameter read in from input csv or xls
                %strutR=StrutDia/2; % Current strut radius
                %Associated Node Coordiantes
                X1=Nodes(curStrut(1),1:3); % Point 1 which is the centre of the first node in the connectivity data
                X2=Nodes(curStrut(2),1:3);% Point 2 which is the centre of the second node in the connectivity data
                %% Calculate Strut properties
                strutL=pdist([X1;X2],'euclidean'); % Strut Axial Length calculated via the euclidean between sphere centre positions
                strutAspectRatio=strutL/StrutDia;
                %Angle to Platen
                V1=X2-X1; % Vector 1 - Vector from connectivity points (axial vector of current strut)
                X3=X1-[0 X1(2)+1 0]; %Create a coordinate point above the first coordinate X1
                V2=X3-X1; % Vector 2 - Vector from point X1 to X3, this vector is perpendicular to the "platen"
                ang=abs(90-abs(rad2deg(acos(dot(V1, V2) / (norm(V1)*norm(V2)))))); % Calculate smallest angle between vectors 1 & 2 which equals angle to platen
                % Analyse Additive Manufacture Build Quality
                if ang>=40
                    BuildQuality='Robust zone';
                    RobustZoneCounter=RobustZoneCounter+1;
                elseif (30<=ang)&&(ang<=40)
                    BuildQuality='Compromised zone';
                    CompromisedZoneCounter=CompromisedZoneCounter+1;
                elseif ang<=30
                    BuildQuality='Failed Zone';
                    FailedZoneCounter=FailedZoneCounter+1;
                    
                end
                % Write current strut details to summaryData cell array for all strut
                % details for excel spreadsheet output
                summaryData{idxstrut,1}=idxstrut;
                summaryData{idxstrut,2}=curStrut(1);
                summaryData{idxstrut,3}=curStrut(2);
                summaryData{idxstrut,4}=curStrut(3);
                summaryData{idxstrut,5}=strutL;
                summaryData{idxstrut,6}=ang;
                summaryData{idxstrut,7}=BuildQuality;
                summaryData{idxstrut,8}=strutAspectRatio;
                
            end
            
            %% plot angles
            angles = cell2mat(summaryData(:,6));
            numColours = 18;
            colours = interp1(linspace(0,90,numColours),colormap(jet(numColours)),angles,'nearest');
            plot(obj,colours);
            a = gca;
            a.CLim = [0,90];
            colormap(jet(numColours))
            colorbar;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            %% Calculate Redundancy
            strutTest=[cell2mat(summaryData(:,6)) cell2mat(summaryData(:,8))];
            RowData=uniquetol(strutTest(:,1),0.05);
            summaryRedun=[];
            for idxRed=1:length(RowData)
                curAngle=RowData(idxRed);
                [cnt,~]=find(strutTest(:,1)<curAngle+0.05 & strutTest(:,1)>curAngle-0.05);
                summaryRedun=[summaryRedun; curAngle length(cnt)];
            end
            %% CREATE STRUT PROPERTIES REVIEW XLS
            Headers={'Maxwell Num','Num.of Robust Zone','Num.of Compromised Zone','Num.of Failed Zone';LatBehavOut,RobustZoneCounter,CompromisedZoneCounter,FailedZoneCounter};
            Headers2={'Strut','Node 1','Node 2','Strut Diameter','Strut Length','Angle to Platen','Strut Quality Manufacture Zone','Aspect Ratio'};
            xlswrite(fileName,Headers);
            xlswrite(fileName,Headers2,'Sheet1','A3');
            xlswrite(fileName,summaryData,'Sheet1','A4');
            % Write repeated struts
            xlswrite(fileName,{'Repeated Strut Details'},'Sheet1','J2');
            xlswrite(fileName,{'Strut Angle'},'Sheet1','J3');
            xlswrite(fileName,{'Repetitions'},'Sheet1','K3');
            xlswrite(fileName,summaryRedun,'Sheet1','J4');
        end
        function obj = load(obj,file)
            % load Loads data for use in PLG
            %    file is passed from main PLG.
            %    Returns a modified PLG object.
            parts = strsplit(file,'.');
            extension = parts{end};
            switch extension
                case obj.loadExtensions{1,1}
                    % xml - use xml2lattice first or use as a unit cell
                    error('PLG:load','XML is not a suitable load format use xml2lattice if desired.');
                case obj.loadExtensions{2,1}
                    % stl - use stlHandler class
                    error('PLG:load','A stl file is not a suitable input for the PLG function use stlHandler.');
                case obj.loadExtensions{3,1}
                    % lattice - assumes lattice model type
                    data = csvread(file);
                case obj.loadExtensions{4,1}
                    % csv - assumes lattice model type
                    warning('PLG:load','A csv will be loaded assuming it is the same format as a .lattice file.');
                    data = csvread(file);
                case obj.loadExtensions{5,1}
                    % excel - assumes lattice model
                    warning('PLG:load','A excel will be loaded assuming it is the same format as a .lattice file.');
                    data = xlsread(file);
                otherwise
                    error('PLG:load','Not a suitable load format, a .lattice file is prefered.');
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
    end          
    methods % save out methods
        function saveStl(obj,fullName)
            % saveStl Save out an stl file
            
            if isempty(obj.resolution) || isempty(obj.sphereResolution)
                error('please set resolution');
            end
            
            numFacets = obj.resolution*4;%number of facets created for one strut
            numLinks = size(obj.struts,1);
            numVertices = size(obj.vertices,1);
            totalFacetsNoBall = numFacets*numLinks;
            if all(obj.sphereDiameter==0) || ~obj.sphereAddition
                % do not write spheres
                totalFacets = totalFacetsNoBall;
            else
                numSpheres = sum(obj.sphereDiameter~=0);
                totalFacets = totalFacetsNoBall + 2*obj.sphereResolution*(obj.sphereResolution-1)*numSpheres;
            end
            %write the header
            fid=fopen(fullName,'w');
            fprintf(fid, '%-80s', 'fast stl generator'); %binary write information
            fwrite(fid,totalFacets,'uint32'); %stl binary header file contains the total number of facets in the stl file
            
            % write the facets
            %write out the struts
            for inc = 1:numLinks
                writeSingleStrut(obj,fid,inc); % adds a single strut to the stl file
            end
            
            if all(obj.sphereDiameter==0) || ~obj.sphereAddition
                fclose(fid);
                return;
            end
            %write the spheres
            [x,y,z]=sphere(obj.sphereResolution); %create sphere 
            ball.struts= convhull([x(:), y(:), z(:)]); %create triangle links
            sizer = size(ball.struts,1);
            ball.vertices=[x(:),y(:),z(:)]; %store the points
            if obj.baseFlat
                shouldBeFlat = abs(obj.vertices(:,3)-min(obj.vertices(:,3)))<1e-5;
                test = z<0;
                z(test)=-1;
                flatBase.struts = ball.struts;
                flatBase.vertices = [x(:),y(:),z(:)];
            else
                shouldBeFlat=false(numVertices,1);
            end
            
            for inc = 1:numVertices
                if obj.sphereDiameter(inc)~=0
                    if ~shouldBeFlat(inc)
                        writeSingleSphere(obj,fid,ball,sizer,inc); % adds a single sphere to the stl file
                    else
                        writeSingleSphere(obj,fid,flatBase,sizer,inc); % adds a flat based sphere to the stl file
                    end
                    
                end
            end
            fclose(fid);
        end
        function saveAbaqus(obj,fullName)
            % saveAbaqus Save out an abaqus beam model.
            
            fid = fopen(fullName,'w');
            
            % Write the header
            fprintf(fid,'*Heading\n');
            fprintf(fid,'** Job name: Job-1 Model name: Model-1\n');
            fprintf(fid,'** Generated by: Programatic Lattice generator\n');
            fprintf(fid,'*Preprint, echo=NO, model=NO, history=NO, contact=NO\n');
            fprintf(fid,'** PARTS\n');
            fprintf(fid,'*Part, name=Lattice_%s_%dx%dx%d\n',obj.latticeType,obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
            
            % write the nodes
            fprintf(fid,'*Node\n');
            numNodes=size(obj.vertices,1);
            for inc = 1:numNodes
                x = obj.vertices(inc,1);
                y = obj.vertices(inc,2);
                z = obj.vertices(inc,3);
                fprintf(fid,'\t%d,\t%0.6e,\t%0.6e,\t%0.6e\n',inc,x,y,z);
            end
            
            % write the connections
            fprintf(fid,'*Element, type=B31\n');
            numLinks=size(obj.struts,1);
            for inc = 1:numLinks
                c1 = obj.struts(inc,1);
                c2 = obj.struts(inc,2);
                fprintf(fid,'\t%d,\t%0.6e,\t%0.6e\n',inc,c1,c2);
            end
            
            fclose(fid);
        end
        function saveLattice(obj,fullName)
            % saveLattice saves data to a csv file with the .lattice 
            % extension
            %    fullname is the path.
            %    Returns void.
            
            % write out the data
            numNodes=size(obj.vertices,1);
            numLinks=size(obj.struts,1);
            
            dlmwrite(fullName,numNodes);
            dlmwrite(fullName,numLinks,'-append');
            data = [obj.vertices,obj.sphereDiameter];
            dlmwrite(fullName,data,'-append');
            
            data = [obj.struts,obj.strutDiameter];
            dlmwrite(fullName,data,'-append');
        end
        function save3mf(obj,fullName)
            % save3mf saves data to a 3mf format
            %    fullname is the path.
            %    Returns void.
            %
            % Save the PLG lattice as the 3D manufacturing file format
            % This consists of three sections
            % write a xml with the following general structure
            %  model
            %    resources
            %      object
            %        mesh
            %           vertices
            %             vertex
            %           triangles
            %             triangle
            %        components
            %           component
            %   build
            %
            % 1. write a generic strut
            % 2. write a generic sphere
            % 3. write all struts and spheres using a series of transformations
            
            % setup object ID
            ID.strut = '0'; % generic strut
            ID.ball =  '1';
            if obj.baseFlat
                ID.flatBall = '2';
                ID.componentsStrut = '3';
                ID.componentsBall = '4';
                ID.componentsFlatBall = '5';
            else
                ID.componentsStrut = '2';
                ID.componentsBall = '3';
            end
            
            %create an object to hold xml and then generate its structure
            xmlObject = com.mathworks.xml.XMLUtils.createDocument('model');
            model = xmlObject.getDocumentElement;
            model.setAttribute('unit','millimeter');
            model.setAttribute('xml:lang','en-US');
            
            resources = xmlObject.createElement('resources');
            model.appendChild(resources);
            
            strutObject = xmlObject.createElement('object');
            strutObject.setAttribute('id',ID.strut);
            strutObject.setAttribute('type','model');
            resources.appendChild(strutObject);
            
            ballObject = xmlObject.createElement('object');
            ballObject.setAttribute('id',ID.ball);
            ballObject.setAttribute('type','model');
            resources.appendChild(ballObject);
            
            strutMesh = xmlObject.createElement('mesh');
            strutObject.appendChild(strutMesh);
            
            ballMesh = xmlObject.createElement('mesh');
            ballObject.appendChild(ballMesh);
            
            strutVertices = xmlObject.createElement('vertices');
            strutMesh.appendChild(strutVertices);
            strutTriangles = xmlObject.createElement('triangles');
            strutMesh.appendChild(strutTriangles);
            
            ballVertices = xmlObject.createElement('vertices');
            ballMesh.appendChild(ballVertices);
            ballTriangles = xmlObject.createElement('triangles');
            ballMesh.appendChild(ballTriangles);
            
            % write the individual strut and ball
            xmlObject = threeMfStrut(obj,xmlObject,strutVertices,strutTriangles);
            if ~obj.sphereAddition
                obj.sphereResolution = 4;
            end
            xmlObject = threeMfBall(obj,xmlObject,ballVertices,ballTriangles);
            % if required generate the flat base ball and renumber
            if obj.baseFlat
                % use a flat base
                flatBallObject = xmlObject.createElement('object');
                flatBallObject.setAttribute('id',ID.flatBall);
                flatBallObject.setAttribute('type','model');
                resources.appendChild(flatBallObject);
                
                flatBallMesh = xmlObject.createElement('mesh');
                flatBallObject.appendChild(flatBallMesh);
                
                flatBallVertices = xmlObject.createElement('vertices');
                flatBallMesh.appendChild(flatBallVertices);
                flatBallTriangles = xmlObject.createElement('triangles');
                flatBallMesh.appendChild(flatBallTriangles);
                
                % write the flat ball object
                xmlObject = threemfFlatBall(obj,xmlObject,flatBallVertices,flatBallTriangles);
            end

            
            % write out the replications that goes under components
            %replicate the struts throughout space
            strutObject1 = xmlObject.createElement('object');
            strutObject1.setAttribute('id',ID.componentsStrut);
            strutObject1.setAttribute('type','model');
            resources.appendChild(strutObject1);
            strutComponent = xmlObject.createElement('components');
            strutObject1.appendChild(strutComponent);
            
            xmlObject = threeMfStrutReplicate(obj,xmlObject,strutComponent,ID.strut);
            
            if obj.sphereAddition
                ballObject1 = xmlObject.createElement('object');
                ballObject1.setAttribute('id',ID.componentsBall);
                ballObject1.setAttribute('type','model');
                resources.appendChild(ballObject1);
                ballComponent = xmlObject.createElement('components');
                ballObject1.appendChild(ballComponent);
                if obj.baseFlat
                    flatBallObject1 = xmlObject.createElement('object');
                    flatBallObject1.setAttribute('id',ID.componentsFlatBall);
                    flatBallObject1.setAttribute('type','model');
                    resources.appendChild(flatBallObject1);
                    flatBallComponent = xmlObject.createElement('components');
                    flatBallObject1.appendChild(flatBallComponent);
                    
                    xmlObject = threeMfFlatBallReplicate(obj,xmlObject,ballComponent,flatBallComponent,ID.ball,ID.flatBall);
                else
                    xmlObject = threeMfBallReplicate(obj,xmlObject,ballComponent,ID.ball);
                end
            end
            
            % write out the build section which states that the components
            % should be put in the assembly
            build = xmlObject.createElement('build');
            model.appendChild(build);
            item = xmlObject.createElement('item');
            item.setAttribute('objectid',ID.componentsStrut);
            build.appendChild(item);
            if obj.sphereAddition
                item = xmlObject.createElement('item');
                item.setAttribute('objectid',ID.componentsBall);
                build.appendChild(item);
                if obj.baseFlat
                    item = xmlObject.createElement('item');
                    item.setAttribute('objectid',ID.componentsFlatBall);
                    build.appendChild(item);
                end
            end
            % gather the other supplementary files that are required to make a 3mf file
            mkdir('_rels');
            mkdir('3D');
            xmlwrite(['3D',filesep,'3dmodel.model'],xmlObject);
            copyfile([obj.runLocation,filesep,'for_3mf/.rels'],'_rels/.rels');
            copyfile([obj.runLocation,filesep,'for_3mf/[Content_Types].xml'],'[Content_Types].xml');
            zip('out',{'_rels','3D','[Content_Types].xml'});
            movefile('out.zip',fullName,'f');
            
            % delete the unzipped files
            delete('[Content_Types].xml');
            rmdir('3D','s');
            rmdir('_rels','s');
        end
    end
    methods (Access=protected) % stl format
        function writeSingleStrut(obj,fid,inc)
            % based on a strut increment write a single strut to file 
            point1 = obj.vertices(obj.struts(inc,1),:);
            point2 = obj.vertices(obj.struts(inc,2),:);
            
            vector = point2-point1;
            u1 = vector/norm(vector);
            if u1(3)==1 || u1(3)==-1
                v1 = [1,0,0];
            else
                v1 = cross([0,0,1],u1);
                v1 = v1/norm(v1);
            end
            
            % generate all the points on the strut
            offset = obj.strutDiameter(inc)*v1/2; % has to use radius
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
            
            % write the points to file in the write order to create a stl file
            for j=1:obj.resolution
                % end of strut at point 1
                datOut = circshift(vert1end,j);
                facet_a=point1; facet_b=datOut(1,:); facet_c=datOut(2,:);
                writeSingleFace(obj,fid,facet_a,facet_c,facet_b);
                
                % end of strut at point 2
                datOut = circshift(vert2end,j);
                facet_a=point2; facet_b=datOut(1,:); facet_c=datOut(2,:);
                writeSingleFace(obj,fid,facet_a,facet_c,facet_b);
                
                % along direction point 1 to point 2
                datOut1 = circshift(vert1end,j);
                datOut2 = circshift(vert2end,j);
                facet_a=datOut1(1,:); facet_b=datOut2(1,:); facet_c=datOut2(2,:);
                writeSingleFace(obj,fid,facet_a,facet_c,facet_b);
                
                % along direction point 2 to point 1
                datOut1 = circshift(vert1end,j);
                datOut2 = circshift(vert2end,j);
                facet_a=datOut1(1,:); facet_b=datOut1(2,:); facet_c=datOut2(2,:);
                writeSingleFace(obj,fid,facet_a,facet_b,facet_c);
            end
        end
        function writeSingleSphere(obj,fid,ball,sizer,inc)
            offset=ball.vertices*obj.sphereDiameter(inc)/2;
            target=[offset(:,1)+obj.vertices(inc,1),offset(:,2)+obj.vertices(inc,2),offset(:,3)+obj.vertices(inc,3)];
            for j=1:sizer
                facet_a=target(ball.struts(j,1),:); facet_b=target(ball.struts(j,2),:); facet_c=target(ball.struts(j,3),:);
                writeSingleFace(obj,fid,facet_a,facet_b,facet_c);
            end
        end
        function writeSingleFace(obj,fid,v1,v2,v3)
            normal=cross(v2-v1,v3-v1);
            fwrite(fid,normal,'float32');           % write normal vector floating point numbers
            fwrite(fid,v1,'float32');   % first vertex (x,y,z)
            fwrite(fid,v2,'float32');   % second vertex
            fwrite(fid,v3,'float32');   % Third vertex
            fwrite(fid,0,'uint16','l');
        end
    end
    methods (Access=protected)% 3mf format
        function xmlObject = threeMfBall(obj,xmlObject,ballVertices,ballTriangles)
            % write a single ball as an object
            % write the vertexs
            [x,y,z]=sphere(obj.sphereResolution); %create sphere with higher accuracy
            ball.struts= convhull([x(:), y(:), z(:)]); %create triangle links
            ball.vertices=[x(:),y(:),z(:)]; %store the points
            ball.vertices = ball.vertices;
            for inc = 1:size(ball.vertices,1)
                vertexNode = xmlObject.createElement('vertex');
              
                val = sprintf('%3.8f',ball.vertices(inc,1));
                vertexNode.setAttribute('x',val);
                val = sprintf('%3.8f',ball.vertices(inc,2));
                vertexNode.setAttribute('y',val);
                val = sprintf('%3.8f',ball.vertices(inc,3));
                vertexNode.setAttribute('z',val);
                ballVertices.appendChild(vertexNode);
            end
            
            %write the connections
            for inc = 1:size(ball.struts,1)
                triangleNode = xmlObject.createElement('triangle');
                
                val = sprintf('%0.0f',ball.struts(inc,1)-1);
                triangleNode.setAttribute('v1',val);
                val = sprintf('%0.0f',ball.struts(inc,2)-1);
                triangleNode.setAttribute('v2',val);
                val = sprintf('%0.0f',ball.struts(inc,3)-1);
                triangleNode.setAttribute('v3',val);
                
                ballTriangles.appendChild(triangleNode);
            end
        end
        function xmlObject = threemfFlatBall(obj,xmlObject,ballVertices,ballTriangles)
            % write a single ball as an object
            % write the vertexs
            [x,y,z]=sphere(obj.sphereResolution); %create sphere with higher accuracy
            test = z<0;
            z(test)=-1;
            ball.struts= convhull([x(:), y(:), z(:)]); %create triangle links
            ball.vertices=[x(:),y(:),z(:)]; %store the points
            ball.vertices = ball.vertices;
            for inc = 1:size(ball.vertices,1)
                vertexNode = xmlObject.createElement('vertex');
              
                val = sprintf('%3.8f',ball.vertices(inc,1));
                vertexNode.setAttribute('x',val);
                val = sprintf('%3.8f',ball.vertices(inc,2));
                vertexNode.setAttribute('y',val);
                val = sprintf('%3.8f',ball.vertices(inc,3));
                vertexNode.setAttribute('z',val);
                ballVertices.appendChild(vertexNode);
            end
            
            %write the connections
            for inc = 1:size(ball.struts,1)
                triangleNode = xmlObject.createElement('triangle');
                
                val = sprintf('%0.0f',ball.struts(inc,1)-1);
                triangleNode.setAttribute('v1',val);
                val = sprintf('%0.0f',ball.struts(inc,2)-1);
                triangleNode.setAttribute('v2',val);
                val = sprintf('%0.0f',ball.struts(inc,3)-1);
                triangleNode.setAttribute('v3',val);
                
                ballTriangles.appendChild(triangleNode);
            end
            
        end
        function xmlObject = threeMfStrut(obj,xmlObject,strutVertices,strutTriangles)
            % write a single strut as an object
            radius = 0.5; % as the transform contains the scaling for the strut diameter
            point1 = [0,0,-0.5];
            point2 = [0,0,0.5];
            
            vector = point2-point1;
            u1 = vector/norm(vector);
            crosser = [0,1,0];
            v1 = cross(crosser,u1);
            v1 = v1/norm(v1);
            offset = radius*v1;
                
            % verts
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
            vertsOut = [point1;vert1end;point2;vert2end];
                
            %struts
            s1 = zeros(obj.resolution,3);
            s2 = zeros(obj.resolution,3);
            s3 = zeros(obj.resolution,3);
            s4 = zeros(obj.resolution,3);
            end1 = [1:obj.resolution]+1;
            end2 = [1:obj.resolution]+obj.resolution+2;
            end1 = circshift(end1,-1,2);
            end2 = circshift(end2,-1,2);
            for j=1:obj.resolution
                end1 = circshift(end1,1,2);
                end2 = circshift(end2,1,2);
                %middle point1 -> end 1
                s1(j,:) = [1,end1(2),end1(1)];
                %middle point2 -> end 2
                s2(j,:) = [obj.resolution+2,end2(1),end2(2)];
                %end1 -> end2
                s3(j,:) = [end1(1),end1(2),end2(1)];
                %end2 -> end1
                s4(j,:) = [end2(2),end2(1),end1(2)];
            end
            strutCons = [s1;s2;s3;s4];
            
            % write verts then connections
            for inc = 1:size(vertsOut,1)
                vertexNode = xmlObject.createElement('vertex');
                
                val = sprintf('%3.8f',vertsOut(inc,1));
                vertexNode.setAttribute('x',val);
                val = sprintf('%3.8f',vertsOut(inc,2));
                vertexNode.setAttribute('y',val);
                val = sprintf('%3.8f',vertsOut(inc,3));
                vertexNode.setAttribute('z',val);
                
                strutVertices.appendChild(vertexNode);
            end

            for inc = 1:size(strutCons,1)
                triangleNode = xmlObject.createElement('triangle');
                
                val = sprintf('%0.0f',strutCons(inc,1)-1);
                triangleNode.setAttribute('v1',val);
                val = sprintf('%0.0f',strutCons(inc,2)-1);
                triangleNode.setAttribute('v2',val);
                val = sprintf('%0.0f',strutCons(inc,3)-1);
                triangleNode.setAttribute('v3',val);
                
                strutTriangles.appendChild(triangleNode);
            end
        end
        function xmlObject = threeMfStrutReplicate(obj,xmlObject,strutComponent,copyObjId)
            % for every strut generate a transform and then write out said
            % transform
            numStruts = size(obj.struts,1);
            for inc = 1:numStruts
                % convert a strut to a transform
                currentTransform = beam2tramsform(obj,inc);
                str = sprintf('%5.5f ',currentTransform);
                str(end) = [];
                
                componentNode = xmlObject.createElement('component');
                componentNode.setAttribute('objectid',copyObjId);
                componentNode.setAttribute('transform',str);
                strutComponent.appendChild(componentNode);
            end
        end
        function xmlObject = threeMfBallReplicate(obj,xmlObject,ballComponent,copyObjId)
            % write the ball locations
            for inc = 1:length(obj.sphereDiameter)
                point = obj.vertices(inc,:);
                currentTransform = eye(4)*obj.sphereDiameter(inc)/2;
                currentTransform(4,1:3) = point;
                currentTransform(:,4) = [];
                str = sprintf('%5.5f ',currentTransform');
                str(end) = [];
                
                componentNode = xmlObject.createElement('component');
                componentNode.setAttribute('objectid',copyObjId);
                componentNode.setAttribute('transform',str);
                ballComponent.appendChild(componentNode);
            end
        end
        function xmlObject = threeMfFlatBallReplicate(obj,xmlObject,ballComponent,flatBallComponent,copyObjIdBall,copyObjIdFlat)
            % write the ball locations
            minZ = min(obj.vertices(:,3));
            for inc = 1:length(obj.sphereDiameter)
                point = obj.vertices(inc,:);
                currentTransform = eye(4)*obj.sphereDiameter(inc)/2;
                currentTransform(4,1:3) = point;
                currentTransform(:,4) = [];
                str = sprintf('%5.5f ',currentTransform');
                str(end) = [];
                if abs(minZ-point(3))<1e-5
                    % use flat ball
                    componentNode = xmlObject.createElement('component');
                    componentNode.setAttribute('objectid',copyObjIdFlat);
                    componentNode.setAttribute('transform',str);
                    flatBallComponent.appendChild(componentNode);
                else
                    % use normal ball
                    componentNode = xmlObject.createElement('component');
                    componentNode.setAttribute('objectid',copyObjIdBall);
                    componentNode.setAttribute('transform',str);
                    ballComponent.appendChild(componentNode);
                end
            end
        end
        function transform = beam2tramsform(obj,inc)
            % beam2tramsform
            %    Determine the affine transform for a given beam increment.
            verts = [0,0,-0.5;0,0,0.5];
            thirdPoint = verts(1,:)+[0.5,0,0];
            forthPoint = verts(2,:)+[0,0.5,0];
            originalPoints = [verts(1,:),1;verts(2,:),1;thirdPoint,1;forthPoint,1];
            
            point1 = obj.vertices(obj.struts(inc,1),:);
            point2 = obj.vertices(obj.struts(inc,2),:);
            vector = point2-point1;
            u = vector/norm(vector);
            if abs(u(3))==1
                    crosser = [1,0,0];
                else % perfect z or any other vector
                    crosser = [0,0,1];
            end
            v = cross(crosser,u);
            v = v/norm(v);
            offset = obj.strutDiameter(inc)*v/2;
            
            point3 = point1+offset;
            vector = point3-point1;
            w = vector/norm(vector);
            x = cross(w,u);
            x = x/norm(x);
            offset = obj.strutDiameter(inc)*x/2;
            point4 = point2-offset;
            
            newPoints = [point1,1;point2,1;point3,1;point4,1;];
            affine = originalPoints\newPoints;
            transform = [affine(1,1:3),affine(2,1:3),affine(3,1:3),affine(4,1:3)];
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
