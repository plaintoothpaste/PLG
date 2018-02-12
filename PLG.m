classdef PLG
    % The PLG program rewritten in class form
    % only core lattice generation functions will be stored here
    % an input is required to use this
    % EXAMPLE
    % obj = PLG('bcc',12,0.3,...
    %     1,0.5,8,2,2,2,...
    %     4,4,5,0,0,0);
    % save(obj);
    properties (SetAccess=protected)
        unitName;
        unitType;
        
        resolution;
        strutDiameter;
        
        sphereAddition;
        sphereResolution;
        sphereDiameter;
        
        vertices;
        struts;
        transform; % holds the transform for a specific strut to generate the entire lattice structure
        
        % 3 value vectors for the x y z params respectivally
        unitSize
        replications
        origin
        
        validExtensions  = {'xlsx'; 'csv'; 'custom'};
        filterSpecOut = {'*.stl','3D geometry file';...
            '*.amf','additive manufacturing format';...
            '*.inp','Abaqus input file';...
            '*.bin','Binary storage file';...
            '*.xlsx','Excel format';...
            '*.custom','Custom csv format for debug etc'};
        strutureType;
        tolerance; % defined as 1/100 of the shortest length present
        dx; % length of a strut
    end
    methods
        function obj = PLG(varargin)
            % creates the PLG object
            switch numel(varargin)
                case 0
                    % generate a new lattice
                    obj.sphereAddition = false;
                case 1
                    % import a custom lattice file containing beam and node
                    % definitions see load function for more information
                    obj = load(obj,varargin{1});
                    obj.strutureType = 0;
                case 15
                    %Generate a new regular lattice batch methods (legacy)
                    %   input order same as properties order
                    obj.latticeType = varargin{1};
                    obj.resolution = varargin{2};
                    obj.strutDiameter = varargin{3};
                    if varargin{4}==1
                        obj.sphereDiameter = varargin{5};
                        obj.sphereResolution = varargin{6};
                    else
                        obj.sphereDiameter = [];
                        obj.sphereResolution = 2;
                    end
                    
                    obj.unitSize(1) = varargin{7};
                    obj.unitSize(2) = varargin{8};
                    obj.unitSize(3) = varargin{9};
                    obj.replications(1) = varargin{10};
                    obj.replications(2) = varargin{11};
                    obj.replications(3) = varargin{12};
                    obj.origin(1) = varargin{13};
                    obj.origin(2) = varargin{14};
                    obj.origin(3) = varargin{15};
                    
                    % set tolerance as required in clean lattice
                    obj.tolerance = min([obj.unitSize(1),obj.unitSize(2),obj.unitSize(3)])/100;
                    
                    % generate the lattice structure based on inputs
                    obj = latticeGenerate(obj); % generate structure
                    
                otherwise
                    error('Incorrect number of inputs');
            end
        end
        function obj = set(obj,varargin)
            % set a value in the PLG that can be edited
            allowable={'resolution','strutDiameter','sphereAddition','sphereResolution','sphereDiameter',...
                'unitSize','replications','origin'};
            dataType = {'double','double','logical','double','double','double','double','double'};
            dataSize = {[1,1],[1,1],[1,1],[1,1],[1,1],[1,3],[1,3],[1,3]};
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
        function obj = defineUnit(obj,unitNames,type)
            % define a unit cell and wheter output is radial or cartesian
            addpath('unitCell');
            if ~exist('type','var')
                type = 'cartesian';
                warning('Coordinate system not specified using cartesian')
            end
            if ~iscell(unitNames)
                error('unitNames input must be a cell and can contain multiple names')
            end
            
            
            switch type
                case 'cartesian'
                    % for the list of unit cell names generate an object that holds them all
                    if length(unitNames)>1
                        unitObj = unitCell(unitNames{1});
                        for inc = 2:length(unitNames)
                            unitObj = addUnit(unitObj,unitNames{inc});
                        end
                    else
                        unitObj = unitCell(unitNames{1});
                    end
                    unitObj = set(unitObj,'diameter',obj.strutDiameter,'resolution',obj.resolution,'scale',obj.unitSize);
                    [obj.vertices,obj.struts,obj.transform,obj.unitName,obj.unitType] = output(unitObj);
                case 'radial'
                    error('radial/polar coordinate system not yet supported');
                otherwise
                    error('coordinate system type not supported');
                    
            end
            rmpath('unitCell');
        end
                
        function obj = cellReplication(obj)
            % if unitType is a beam it will replicate transformation if it is a facet type it will
            % replicate struts and verts
            % replaces obj.replications with a list of x,y,z coordinates for the unit cell
            xPlacement = 0:obj.unitSize(1):obj.unitSize(1)*(obj.replications(1)-1);
            yPlacement = 0:obj.unitSize(2):obj.unitSize(2)*(obj.replications(2)-1);
            zPlacement = 0:obj.unitSize(3):obj.unitSize(3)*(obj.replications(3)-1);
            [XX,YY,ZZ] = ndgrid(xPlacement,yPlacement,zPlacement);
            
            switch obj.unitType
                case 'beam'
                    % replicate transform
                    originalTransforms = obj.transform;
                    sizeTransforms = size(originalTransforms);
                    numReps = length(XX(:));
                    obj.transform = zeros(numReps*sizeTransforms(1),sizeTransforms(2));
                    for incReps = 1:numReps
                        newX = XX(incReps);
                        newY = YY(incReps);
                        newZ = ZZ(incReps);
                        for incTrans = 1:sizeTransforms(1)
                            currentInc = incTrans + (incReps-1)*sizeTransforms(1);
                            currentTransform = originalTransforms(incTrans,:);
                            affine = [currentTransform(1:3),0;currentTransform(4:6),0;currentTransform(7:9),0;currentTransform(10:12),1];
                            newAffine = zeros(4); newAffine(4,1:3) = [newX,newY,newZ];
                            newAffine = newAffine+affine;
                            obj.transform(currentInc,:) = [newAffine(1,1:3),newAffine(2,1:3),newAffine(3,1:3),newAffine(4,1:3)];
                        end
                    end
                case 'facet'
                    vertOut = arrayfun(@(x,y,z) ...
                        [obj.vertices(:,1) + x,...
                        obj.vertices(:,2) + y,...
                        obj.vertices(:,3) + z]...
                        ,XX,YY,ZZ,'UniformOutput',0);
                    obj.vertices = cell2mat(vertOut(:));
                    
                    strutCounter = 0:(numel(XX)-1);
                    numStruts = max(max(obj.struts));
                    strutCounter = strutCounter*numStruts;
                    strutOut = arrayfun(@(x) obj.struts + x,strutCounter,'UniformOutput',0);
                    obj.struts = cell2mat(strutOut(:));
                otherwise
                    error('cell replication can not be performed until a unit cell is defined or multiple times');
            end
        end
        function obj = cleanLattice(obj)
            % cleans the lattice structure including the removal of duplicate transforms or
            % duplicate vertices and struts depending on type
            switch obj.unitType
                case 'beam'
                    % determine identical transforms and remove can be done with a unique
                    tester = round(obj.transform,6,'significant');
                    [~,I] = unique(tester,'rows');
                    obj.transform = obj.transform(I,:);
                case 'facet'
                    % tolerance - get the shortest strut and divide by 1000
                    verts1 = obj.vertices(obj.struts(:,1),:);
                    verts2 = obj.vertices(obj.struts(:,2),:);
                    lengthVerts = sum(sqrt((verts1-verts2).^2),2);
                    obj.tolerance = min(lengthVerts)/15;
                    % duplicate vertices
                    [obj.vertices,i,indexn]=uniquetol(obj.vertices,obj.tolerance,'ByRows',1,'DataScale',1);
                    if ~isempty(obj.sphereDiameter)
                        obj.sphereDiameter = obj.sphereDiameter(i);
                    end
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
                otherwise
                    error('cleaning  can not be performed on this unitType');
            end
        end
        
        function obj = custom2beam(obj)
            % changes a custom imported lattice file to a series of
            % transforms this enables 3mf and amf save out as well as more
            % efficient storage.
            
            error('TODO');
        end
        function plot(obj,colours)
            % plot the lattice with nodes highlighted currently only works for beams
            switch obj.unitType
                case 'beam'
                    f = figure;
                    f.Units	= 'Normalized';
                    f.Position = [0,0,1,1];
                    f.Name = 'STL plotter';
                    a = axes;
                    a.View = [45,45];
                    axis vis3d
                    a.NextPlot='add';
                    
                    numTransforms=size(obj.transform,1);
                    if ~exist('colours','var')
                        colours = repmat([0.3,0.3,0.3,0.5],numTransforms,1);
                    end
                    zeroPoints = [obj.vertices,[1;1]];
                    for inc = 1:numTransforms
                        currentTransform = obj.transform(inc,:);
                        affine = zeros(4);
                        affine(4,4) = 1;
                        affine(1,1:3) = currentTransform(1:3);
                        affine(2,1:3) = currentTransform(4:6);
                        affine(3,1:3) = currentTransform(7:9);
                        affine(4,1:3) = currentTransform(10:12);
                        
                        points = zeroPoints*affine;
                        p = plot3(points(:,1),points(:,2),points(:,3),'Color',colours(inc,:));
                        p.MarkerFaceColor = [0.9,0.5,0]; p.MarkerEdgeColor = 'none';
                    end
                case 'facet'
                    
                otherwise
                    error('Plotting not yet supported for this unitType')
            end
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
        
        function obj = translate(obj,x,y,z)
            %translates a lattice in space
            obj.vertices(:,1) = obj.vertices(:,1)+x;
            obj.vertices(:,2) = obj.vertices(:,2)+y;
            obj.vertices(:,3) = obj.vertices(:,3)+z;
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
            
            %rotation x then y then z split for debugging
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
        function obj = plus(obj,obj1)
            % add obj1 to obj2
            newStruts = obj1.struts+max(obj.struts(:));
            obj.struts = [obj.struts;newStruts];
            obj.strutDiameter = [obj.strutDiameter;obj1.strutDiamter];
            obj.sphereDiameter = [obj.sphereDiameter;obj1.sphereDiameter];
            obj.vertices = [obj.vertices;obj1.vertices];
            
            % remove any matching struts
            obj = cleanLattice(obj);
        end
        
        function obj = splitStruts(obj)
            % takes a badly defined input lattice and scoures it for any struts that intersect and
            % but do not connect to each other. It then breaks these struts and creates a new clean
            % lattice
            % this functions should only be used if the user knows what they are doing and may
            % destroy unique diameter data
            maxInd = 0; % current highest index in strutsOut
            lengthstruts = 0; % current length of strutsOut
            numstruts = length(obj.struts);
            vertsOut = [];
            strutsOut = [];
            splitStruts = cell(numstruts,1);
            % get tol
            obj = getTolerance(obj);
            tol = obj.tolerance;
            struts2Check = 1:numstruts;
            boundingBox = PLG.getBound(obj.vertices,obj.struts,struts2Check);
            %% main loop
            for inc=1:numstruts
                % all struts with intersecting BB
                struts2CheckTmp = struts2Check;
                struts2CheckTmp(inc) = [];
                currentBB = boundingBox{inc};
                currentP1 = obj.vertices(obj.struts(inc,1),:);
                currentP2 = obj.vertices(obj.struts(inc,2),:);
                isIntersect = PLG.findBbIntersect(currentBB,boundingBox(struts2CheckTmp));
                potentialStruts = struts2CheckTmp(isIntersect);
                % remove that are already connected properly
                potentialStruts = PLG.removeNormalConStruts(currentP1,currentP2,potentialStruts,obj.struts,obj.vertices,tol);
                potentialStruts = PLG.removeParralelStruts(currentP1,currentP2,potentialStruts,obj.struts,obj.vertices);
                % check if any struts are already split if so grab their sub BB and
                % eleminate any that do not intersect
                % potentialStruts = PLG.findSplitStrutBB(splitStruts,potentialStruts,currentBB,currentP1,currentP2,strutsOut,vertsOut,tol);
                potentialStruts = [potentialStruts',ones(length(potentialStruts),1)];
                if isempty(potentialStruts)
                    % no struts to break up
                    vertsOut = [vertsOut;currentP1;currentP2];
                    strutsOut = [strutsOut;maxInd+1,maxInd+2];
                    maxInd = maxInd+2;
                    lengthstruts = lengthstruts+1; % number of new struts
                else
                    [vertsTmpOut,strutsTmpOut] = PLG.findIntersect(currentP1,currentP2,...
                        potentialStruts,obj.struts,obj.vertices,...
                        strutsOut,vertsOut,tol);
                    vertsOut = [vertsOut;vertsTmpOut];
                    strutsTmpOut = strutsTmpOut+maxInd;
                    % debug: 
                    % patch('struts',obj.struts(inc,:),'Vertices',obj.vertices,'EdgeColor',[1,0,0]);
                    % patch('struts',obj.struts(potentialStruts,:),'Vertices',obj.vertices); hold on;
                    % scatter3(vertsOut(:,1),vertsOut(:,2),vertsOut(:,3));

                    maxInd = strutsTmpOut(end,2);
                    strutsOut = [strutsOut;strutsTmpOut];
                    st = lengthstruts+1; % number of new struts
                    en = st+size(strutsTmpOut,1)-1;
                    splitStruts{inc} = [st,en];
                    lengthstruts = en;
                end
                
            end
            %% final assignment and cleanup note that this destroys strut diameter
            obj.vertices = vertsOut;
            obj.struts = strutsOut;
            % remake sphere and strut diams same as the first value
            obj.strutDiameter = obj.strutDiameter(1);
            obj.sphereDiameter = obj.sphereDiameter(1);
            obj = addDiams(obj);
            obj = cleanLattice(obj);
        end
        function obj = calcDx(obj)
            % returns the absolute vector length
            p1 = obj.vertices(obj.struts(:,1),:);
            p2 = obj.vertices(obj.struts(:,2),:);
            diffx = (p1(:,1)-p2(:,1));
            diffy = (p1(:,2)-p2(:,2));
            diffz = (p1(:,3)-p2(:,3));
            obj.dx = sqrt(diffx.^2+diffy.^2+diffz.^2);
        end
        function getProperties(obj,fileName)
            % convert plg inputs into the function inputs
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
    end
    methods (Access=protected)%not called by the user
        function obj = latticeGenerate(obj)
            % generates a latticeStructure ready for simulation or saving as an stl etc
            % all the below unit cells are stored in their own method section
            switch obj.latticeType
                case 'bcc' % BCC Cell
                    [obj.vertices, obj.struts] = PLG.bcc(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
                case 'bcc2' % BCC Cell with no cantilevers
                    [obj.vertices, obj.struts] = PLG.bcc2(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
                case 'bcz' % BCC Cell
                    [obj.vertices, obj.struts] = PLG.bcz(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
                case 'fcc'
                    [obj.vertices, obj.struts] = PLG.fcc(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
                case 'fccNoXY'
                    [obj.vertices, obj.struts] = PLG.fccNoXY(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
                case 'fbcxyz'
                    [obj.vertices, obj.struts] = PLG.fbcxyz(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
                case 'fcz'
                    [obj.vertices, obj.struts] = PLG.fcz(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
                case 'fbcz'
                    [obj.vertices, obj.struts] = PLG.fbcz(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
                case 'bcc_fcc'
                    [obj.vertices, obj.struts] = PLG.bccFcc(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
                case 'fbc'
                    [obj.vertices, obj.struts] = PLG.fbc(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
                case 'box'
                    [obj.vertices, obj.struts] = PLG.box(obj.unitSize(1),obj.unitSize(2),obj.unitSize(3));
            end
            % translsate the unit cell to the correct location
            obj = translate(obj,obj.origin(1),obj.origin(2),obj.origin(3));
        end
        function obj = load(obj,file)
            % load a custom beam input file to generate a lattice structure
            parts = strsplit(file,'.');
            extension = parts{end};
            switch extension
                case obj.validExtensions{1}
                    % xls
                    data = xlsread(file);
                case obj.validExtensions{2}
                    % csv
                    data = csvread(file);
                case obj.validExtensions{3}
                    % custom
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
    end
    methods % save out methods
        function save(obj)
            % this overloads the save function and allows saving out in various
            % formats however this only saves the latticeStructure structure
            [fileName,pathName,filterIndex] = uiputfile(obj.filterSpecOut);
            file = [pathName,fileName];
            switch filterIndex
                case 1
                    % stl
                    saveStl(obj,file);
                case 2
                    % AMF format
                    saveAmf(obj,file);
                    % TODO
                case 3
                    % ABAQUS
                    saveAbaqus(obj,file);
                case 4
                    % 3mf format
                    save3mf(obj,file)
                case 5
                    saveExcel(obj,file);
                case 6
                    % custom file
                    saveCustom(obj,file)
                otherwise
                    error('No file saved')
            end
        end
        function saveStl(obj,fullName)
            numFacets = obj.resolution*4;%number of facets created for one strut
            numLinks = size(obj.struts,1);
            numVertices = size(obj.vertices,1);
            
            totalFacetsNoBall = numFacets*numLinks;
            totalFacetsWithBall = totalFacetsNoBall + 2*obj.sphereResolution*(obj.sphereResolution-1)*numVertices;
            
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
        function saveAbaqus(obj,fullName)
            % saves out as a abaqus beam model that is used as an input
            % file
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
        function saveExcel(obj,fullName)
            % saves data to an excel compatible format
            numNodes=size(obj.vertices,1);
            numLinks=size(obj.struts,1);
            xlswrite(fullName,numNodes,'Sheet1','A1');
            xlswrite(fullName,numLinks,'Sheet1','A2');
            
            locationVertices=strcat('A3:C',num2str(numNodes+2));
            xlswrite(fullName,obj.vertices,'Sheet1',locationVertices);
            locationstruts=strcat('A',num2str(numNodes+3),':B',num2str(numNodes+2+numLinks));
            xlswrite(fullName,obj.struts,'Sheet1',locationstruts);
            
            % optional write depending on load type
            if obj.strutureType ==0
                % loaded custom
                locationSpheres = strcat('D3:D',num2str(numNodes+2));
                xlswrite(fullName,obj.sphereDiameter,'Sheet1',locationSpheres);
                locationDiameters = strcat('C',num2str(numNodes+3),':C',num2str(numNodes+2+numLinks));
                xlswrite(fullName,obj.strutDiameter,'Sheet1',locationDiameters);
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
                data = ones(numLinks,1)*obj.strutDiameter;
                xlswrite(fullName,data,'Sheet1',locationDiameters);
            end
        end
        function saveCustom(obj,fullName)
            % saves data to an excel compatible format
            numNodes=size(obj.vertices,1);
            numLinks=size(obj.struts,1);
            
            dlmwrite(fullName,numNodes);
            dlmwrite(fullName,numLinks,'-append');
            data = [obj.vertices,obj.sphereDiameter];
            dlmwrite(fullName,data,'-append');
            
            data = [obj.struts,obj.strutDiameter];
            dlmwrite(fullName,data,'-append');
        end
        function saveAmf(obj,fullName)
            % create object to hold the xml data 
            amfDoc = com.mathworks.xml.XMLUtils.createDocument('amf');
            amfNode = amfDoc.getDocumentElement;
            amfNode.setAttribute('unit','millimeter');
            
            %write the geom for a unit cell
            amfDoc = amfUnitObj(obj,amfDoc,amfNode,0);
            
            %write the geom for a single ball
            amfDoc = amfBallObj(obj,amfDoc,amfNode,1);
            
            % replication
            amfDoc = amfReplication(obj,amfDoc,amfNode,0);
            
            xmlwrite(fullName,amfDoc);
            %zip(fullName,fullName);
        end
        function save3mf(obj,fullName)
            % safe a PLG lattice as a 3D manufacturing format file
            if ~strcmp(obj.unitType, 'beam')
                error('only beam unit type can be saved as a 3mf file please use stl or convert to beam');
            end
            
            %create an object to hold xml
            threeMfDoc = com.mathworks.xml.XMLUtils.createDocument('model');
            threeMfNode = threeMfDoc.getDocumentElement;
            threeMfNode.setAttribute('unit','millimeter');
            threeMfNode.setAttribute('xml:lang','en-US');
            
            resourcesNode = threeMfDoc.createElement('resources');
            threeMfNode.appendChild(resourcesNode);
            
            % write the geom for a single strut and ball (if required)
            threeMfDoc = threeMfUnitObj(obj,threeMfDoc,resourcesNode,0);
            if obj.sphereAddition
                threeMfDoc = threeMfBallObj(obj,threeMfDoc,resourcesNode,1);
            end
            %write the replications of the above data
            threeMfDoc = threeMfreplicate(obj,threeMfDoc,threeMfNode,0,1);
            
            % gather the other supplementary files that are required to make a 3mf file
            mkdir('_rels');
            mkdir('3D');
            xmlwrite(['3D',filesep,'3dmodel.model'],threeMfDoc);
            copyfile('other/.rels','_rels/.rels');
            copyfile('other/[Content_Types].xml','[Content_Types].xml');
            zip('out',{'_rels','3D','[Content_Types].xml'});
            movefile('out.zip',[fullName,'.3mf'],'f');
            
            % delete the unzipped files
            delete('[Content_Types].xml');
            rmdir('3D','s');
            rmdir('_rels','s');
        end
    end
    methods % AMF format
        function amfDoc = amfBallObj(obj,amfDoc,amfNode,idNum)
            % write out an object to represent a ball
            objNode = amfDoc.createElement('object');
            objNode.setAttribute('id',num2str(idNum));
            amfNode.appendChild(objNode);
            
            meshNode = amfDoc.createElement('mesh');
            objNode.appendChild(meshNode);
            
            verticesNode = amfDoc.createElement('vertices');
            meshNode.appendChild(verticesNode);
            
            % write the vertexs
            [x,y,z]=sphere(obj.sphereResolution); %create sphere with higher accuracy
            ball.struts= convhull([x(:), y(:), z(:)]); %create triangle links
            ball.vertices=[x(:),y(:),z(:)]; %store the points
            ball.vertices = ball.vertices*obj.strutDiameter/2;
            for inc = 1:size(ball.vertices,1)
                vertexNode = amfDoc.createElement('vertex');
                verticesNode.appendChild(vertexNode);
                
                coordNode = amfDoc.createElement('coordinates');
                vertexNode.appendChild(coordNode);
                
                xNode = amfDoc.createElement('x');
                val = sprintf('%3.8f',ball.vertices(inc,1));
                xNode.appendChild(amfDoc.createTextNode(val));
                coordNode.appendChild(xNode);
                
                yNode = amfDoc.createElement('y');
                val = sprintf('%3.8f',ball.vertices(inc,2));
                yNode.appendChild(amfDoc.createTextNode(val));
                coordNode.appendChild(yNode);
                
                zNode = amfDoc.createElement('z');
                val = sprintf('%3.8f',ball.vertices(inc,3));
                zNode.appendChild(amfDoc.createTextNode(val));
                coordNode.appendChild(zNode);
            end
            
            %write the connections
            volumeNode = amfDoc.createElement('volume');
            meshNode.appendChild(volumeNode);
            for inc = 1:size(ball.struts,1)
                triangleNode = amfDoc.createElement('triangle');
                volumeNode.appendChild(triangleNode);
                
                v1Node = amfDoc.createElement('v1');
                val = sprintf('%0.0f',ball.struts(inc,1)-1);
                v1Node.appendChild(amfDoc.createTextNode(val));
                triangleNode.appendChild(v1Node);
                
                v2Node = amfDoc.createElement('v2');
                val = sprintf('%0.0f',ball.struts(inc,2)-1);
                v2Node.appendChild(amfDoc.createTextNode(val));
                triangleNode.appendChild(v2Node);
                
                v3Node = amfDoc.createElement('v3');
                val = sprintf('%0.0f',ball.struts(inc,3)-1);
                v3Node.appendChild(amfDoc.createTextNode(val));
                triangleNode.appendChild(v3Node);
            end
        end
        function amfDoc = amfUnitObj(obj,amfDoc,amfNode,idNum)
            % write a single unit cell as an object
            radius = obj.strutDiameter/2;
            numLinks = length(obj.struts);
            
            % setup the object
            objNode = amfDoc.createElement('object');
            objNode.setAttribute('id',num2str(idNum));
            amfNode.appendChild(objNode);
            
            meshNode = amfDoc.createElement('mesh');
            objNode.appendChild(meshNode);
            
            verticesNode = amfDoc.createElement('vertices');
            meshNode.appendChild(verticesNode);

            verts = cell(numLinks,1);
            struts = cell(numLinks,1);
            for i=1:numLinks
                point1 = obj.vertices(obj.struts(i,1),:);
                point2 = obj.vertices(obj.struts(i,2),:);
                vector = point2-point1;
                u1 = vector/norm(vector);
                if u1(3)==1 || u1(3)==-1
                    v1 = [1,0,0];
                else
                    v1 = cross([0,0,1],u1);
                    v1 = v1/norm(v1);
                end
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
                verts{i} = [point1;vert1end;point2;vert2end];
                
                %struts
                s1 = zeros(obj.resolution,3);
                s2 = zeros(obj.resolution,3);
                s3 = zeros(obj.resolution,3);
                s4 = zeros(obj.resolution,3);
                end1 = [1:obj.resolution]+1;
                end2 = [1:obj.resolution]+obj.resolution+2;
                end1 = circshift(end1,-1);
                end2 = circshift(end2,-1);
                for j=1:obj.resolution
                    end1 = circshift(end1,1);
                    end2 = circshift(end2,1);
                    %middle point1 -> end 1
                    s1(j,:) = [1,end1(2),end1(1)];
                    %middle point2 -> end 2
                    s2(j,:) = [obj.resolution+2,end2(1),end2(2)];
                    %end1 -> end2
                    s3(j,:) = [end1(1),end1(2),end2(1)];
                    %end2 -> end1
                    s4(j,:) = [end2(2),end2(1),end1(2)];
                end
                struts{i} = [s1;s2;s3;s4];
            end
            
            %% join cells
            strutsOut = [];
            vertsOut = [];
            count = 0;
            for inc = 1:numLinks
                strutsOut = [strutsOut;struts{i}+count];
                count = max(strutsOut(:));
                vertsOut = [vertsOut;verts{inc}];
            end
            
            % make verts unique
            [vertsOut,i,indexn]=uniquetol(vertsOut,1e-8,'ByRows',1,'DataScale',1);
            strutsOut = indexn(strutsOut);
            
            %% write verts then connections
            for inc = 1:size(vertsOut,1)
                vertexNode = amfDoc.createElement('vertex');
                verticesNode.appendChild(vertexNode);
                
                coordNode = amfDoc.createElement('coordinates');
                vertexNode.appendChild(coordNode);
                
                xNode = amfDoc.createElement('x');
                val = sprintf('%3.8f',vertsOut(inc,1));
                xNode.appendChild(amfDoc.createTextNode(val));
                coordNode.appendChild(xNode);
                
                yNode = amfDoc.createElement('y');
                val = sprintf('%3.8f',vertsOut(inc,2));
                yNode.appendChild(amfDoc.createTextNode(val));
                coordNode.appendChild(yNode);
                
                zNode = amfDoc.createElement('z');
                val = sprintf('%3.8f',vertsOut(inc,3));
                zNode.appendChild(amfDoc.createTextNode(val));
                coordNode.appendChild(zNode);
            end
            
            volumeNode = amfDoc.createElement('volume');
            meshNode.appendChild(volumeNode);
            for inc = 1:size(strutsOut,1)
                triangleNode = amfDoc.createElement('triangle');
                volumeNode.appendChild(triangleNode);
                
                v1Node = amfDoc.createElement('v1');
                val = sprintf('%0.0f',strutsOut(inc,1)-1);
                v1Node.appendChild(amfDoc.createTextNode(val));
                triangleNode.appendChild(v1Node);
                
                v2Node = amfDoc.createElement('v2');
                val = sprintf('%0.0f',strutsOut(inc,2)-1);
                v2Node.appendChild(amfDoc.createTextNode(val));
                triangleNode.appendChild(v2Node);
                
                v3Node = amfDoc.createElement('v3');
                val = sprintf('%0.0f',strutsOut(inc,3)-1);
                v3Node.appendChild(amfDoc.createTextNode(val));
                triangleNode.appendChild(v3Node);
            end
        end
        function amfDoc = amfReplication(obj,amfDoc,amfNode,idNum)
            % determine locations for the constellation tool for unit cells
            % then ball
            conNode = amfDoc.createElement('constellation');
            conNode.setAttribute('id',num2str(idNum));
            amfNode.appendChild(conNode);
            
            for incX = 1:obj.replications(1)
                posX = (incX-1)*obj.unitSize(1);
                for incY = 1:obj.replications(2)
                    posY = (incY-1)*obj.unitSize(2);
                    for incZ = 1:obj.replications(3)
                        posZ = (incZ-1)*obj.unitSize(3);
                        
                        instNode = amfDoc.createElement('instance');
                        instNode.setAttribute('objectid','0');
                        conNode.appendChild(instNode);
                        
                        dNode = amfDoc.createElement('deltax');
                        val = sprintf('%0.0f',posX);
                        dNode.appendChild(amfDoc.createTextNode(val));
                        instNode.appendChild(dNode);

                        dNode = amfDoc.createElement('deltay');
                        val = sprintf('%0.0f',posY);
                        dNode.appendChild(amfDoc.createTextNode(val));
                        instNode.appendChild(dNode);
                        
                        dNode = amfDoc.createElement('deltaz');
                        val = sprintf('%0.0f',posZ);
                        dNode.appendChild(amfDoc.createTextNode(val));
                        instNode.appendChild(dNode);
                    end
                end
            end
            
            %balls
            if isempty(obj.sphereDiameter)
                return;
            end
            
            conNode = amfDoc.createElement('constellation');
            conNode.setAttribute('id',num2str(idNum+1));
            amfNode.appendChild(conNode);
            
            obj = cellReplication(obj);
            
            for inc = 1:size(obj.vertices,1)
                instNode = amfDoc.createElement('instance');
                instNode.setAttribute('objectid','1');
                conNode.appendChild(instNode);
                
                dNode = amfDoc.createElement('deltax');
                val = sprintf('%0.0f',obj.vertices(inc,1));
                dNode.appendChild(amfDoc.createTextNode(val));
                instNode.appendChild(dNode);
                
                dNode = amfDoc.createElement('deltay');
                val = sprintf('%0.0f',obj.vertices(inc,2));
                dNode.appendChild(amfDoc.createTextNode(val));
                instNode.appendChild(dNode);
                
                dNode = amfDoc.createElement('deltaz');
                val = sprintf('%0.0f',obj.vertices(inc,3));
                dNode.appendChild(amfDoc.createTextNode(val));
                instNode.appendChild(dNode);
            end
        end
    end
    methods % stl format
        function fid = faceCreate(obj,fid)
            radius = obj.strutDiameter/2;
            numLinks = length(obj.struts);
            
            % calculate points on each facet then write said points to file
            % in the correct order
            for i=1:numLinks
                point1 = obj.vertices(obj.struts(i,1),:);
                point2 = obj.vertices(obj.struts(i,2),:);
                vector = point2-point1;
                u1 = vector/norm(vector);
                if u1(3)==1 || u1(3)==-1
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
                %% join struts
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
            [x,y,z]=sphere(obj.sphereResolution); %create sphere with higher accuracy
            ball.struts= convhull([x(:), y(:), z(:)]); %create triangle links
            sizer = size(ball.struts,1);
            ball.vertices=[x(:),y(:),z(:)]; %store the points
            
            for i=1:size(obj.vertices,1)
                offset=ball.vertices*obj.sphereDiameter(i)/2;
                target=[offset(:,1)+obj.vertices(i,1),offset(:,2)+obj.vertices(i,2),offset(:,3)+obj.vertices(i,3)];
                for j=1:sizer
                    %get values first end
                    facet_a=target(ball.struts(j,1),:);
                    facet_b=target(ball.struts(j,2),:);
                    facet_c=target(ball.struts(j,3),:);
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
    methods % 3mf format
        function threeMfDoc = threeMfBallObj(obj,threeMfDoc,resourcesNode,idNum)
            % write out an object to represent a ball 
            objNode = threeMfDoc.createElement('object');
            objNode.setAttribute('id',num2str(idNum));
            objNode.setAttribute('type','model');
            resourcesNode.appendChild(objNode);
            
            meshNode = threeMfDoc.createElement('mesh');
            objNode.appendChild(meshNode);
            
            verticesNode = threeMfDoc.createElement('vertices');
            meshNode.appendChild(verticesNode);
            
            % write the vertexs
            [x,y,z]=sphere(obj.sphereResolution); %create sphere with higher accuracy
            ball.struts= convhull([x(:), y(:), z(:)]); %create triangle links
            ball.vertices=[x(:),y(:),z(:)]; %store the points
            ball.vertices = ball.vertices*obj.sphereDiameter/2;
            for inc = 1:size(ball.vertices,1)
                vertexNode = threeMfDoc.createElement('vertex');
              
                val = sprintf('%3.8f',ball.vertices(inc,1));
                vertexNode.setAttribute('x',val);
                
                val = sprintf('%3.8f',ball.vertices(inc,2));
                vertexNode.setAttribute('y',val);
                
                val = sprintf('%3.8f',ball.vertices(inc,3));
                vertexNode.setAttribute('z',val);
                
                verticesNode.appendChild(vertexNode);
            end
            
            %write the connections
            triNode = threeMfDoc.createElement('triangles');
            meshNode.appendChild(triNode);
            for inc = 1:size(ball.struts,1)
                triangleNode = threeMfDoc.createElement('triangle');
                
                val = sprintf('%0.0f',ball.struts(inc,1)-1);
                triangleNode.setAttribute('v1',val);
                
                val = sprintf('%0.0f',ball.struts(inc,2)-1);
                triangleNode.setAttribute('v2',val);
                
                val = sprintf('%0.0f',ball.struts(inc,3)-1);
                triangleNode.setAttribute('v3',val);
                
                triNode.appendChild(triangleNode);
            end
        end
        function threeMfDoc = threeMfUnitObj(obj,threeMfDoc,resourcesNode,idNum)
            % write a single unit cell as an object
            radius = obj.strutDiameter/2;
            numLinks = size(obj.struts,1);
            
            % setup the object
            objNode = threeMfDoc.createElement('object');
            objNode.setAttribute('id',num2str(idNum));
            objNode.setAttribute('type','model');
            resourcesNode.appendChild(objNode);
            
            meshNode = threeMfDoc.createElement('mesh');
            objNode.appendChild(meshNode);
            
            verticesNode = threeMfDoc.createElement('vertices');
            meshNode.appendChild(verticesNode);
            

            verts = cell(numLinks,1);
            struts = cell(numLinks,1);
            for i=1:numLinks
                point1 = obj.vertices(obj.struts(i,1),:);
                point2 = obj.vertices(obj.struts(i,2),:);
                vector = point2-point1;
                u1 = vector/norm(vector);
                if u1(3)==1 || u1(3)==-1
                    crosser = [0,1,0];
                    
                else
                    crosser = [1,0,1];
                end
                
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
                verts{i} = [point1;vert1end;point2;vert2end];
                
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
                struts{i} = [s1;s2;s3;s4];
            end
            
            %% join cells
            strutsOut = [];
            vertsOut = [];
            count = 0;
            for inc = 1:numLinks
                strutsOut = [strutsOut;struts{i}+count];
                count = max(strutsOut(:));
                vertsOut = [vertsOut;verts{inc}];
            end
            
            % make verts unique
            [vertsOut,i,indexn]=uniquetol(vertsOut,1e-8,'ByRows',1,'DataScale',1);
            strutsOut = indexn(strutsOut);
            
            %% write verts then connections
            for inc = 1:size(vertsOut,1)
                vertexNode = threeMfDoc.createElement('vertex');
              
                val = sprintf('%3.8f',vertsOut(inc,1));
                vertexNode.setAttribute('x',val);
                
                val = sprintf('%3.8f',vertsOut(inc,2));
                vertexNode.setAttribute('y',val);
                
                val = sprintf('%3.8f',vertsOut(inc,3));
                vertexNode.setAttribute('z',val);
                
                verticesNode.appendChild(vertexNode);
            end
            triNode = threeMfDoc.createElement('triangles');
            meshNode.appendChild(triNode);
            for inc = 1:size(strutsOut,1)
                triangleNode = threeMfDoc.createElement('triangle');
                
                val = sprintf('%0.0f',strutsOut(inc,1)-1);
                triangleNode.setAttribute('v1',val);
                
                val = sprintf('%0.0f',strutsOut(inc,2)-1);
                triangleNode.setAttribute('v2',val);
                
                val = sprintf('%0.0f',strutsOut(inc,3)-1);
                triangleNode.setAttribute('v3',val);
                
                triNode.appendChild(triangleNode);
            end
        end
        function threeMfDoc = threeMfreplicate(obj,threeMfDoc,threeMfNode,idUnit,idBall)
            % for every transform write it to the file
            % build item holds all replications
            buildNode = threeMfDoc.createElement('build');
            threeMfNode.appendChild(buildNode);
            
            numTransform = size(obj.transform,1);
            for inc = 1:numTransform
                itemNode = threeMfDoc.createElement('item');
                itemNode.setAttribute('objectid',num2str(idUnit));
                
                currentTransform = obj.transform(inc,:);
                str = sprintf('%5.5f ',currentTransform);
                str(end) = [];
                itemNode.setAttribute('transform',str);
                buildNode.appendChild(itemNode);
            end
            
            if ~obj.sphereAddition
                return
            end
            % write the ball locations if required
            numTransform =size(obj.transform,1);
            verts = [];
            points = [obj.vertices,[1;1]];
            for inc = 1:numTransform
                currentTransform = obj.transform(inc,:);
                trans = eye(4);
                trans(1,1:3) = currentTransform(1:3);
                trans(2,1:3) = currentTransform(4:6);
                trans(3,1:3) = currentTransform(7:9);
                trans(4,1:3) = currentTransform(10:12);
                newVerts = points*trans;
                verts = [verts;newVerts(:,1:3)];
            end
            verts = unique(verts,'rows');
            for inc = 1:size(verts,1)
                itemNode = threeMfDoc.createElement('item');
                itemNode.setAttribute('objectid',num2str(idBall));
                
                str = sprintf('1 0 0 0 1 0 0 0 1 %5.5f %5.5f %5.5f', verts(inc,1), verts(inc,2),verts(inc,3));
                itemNode.setAttribute('transform',str);
                buildNode.appendChild(itemNode);
            end
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
        function boundingBox = getBound(vertices,struts,desired)
            % get bounding box of ervery strut
            %[minX minY minZ]
            %[maxX maxY maxZ]
            boundingBox = cell(length(desired),1);
            numRequiredBB = length(desired);
            for inc = 1:numRequiredBB
                index = desired(inc);
                maxers = max(vertices(struts(index,1),:),vertices(struts(index,2),:));
                miners = min(vertices(struts(index,1),:),vertices(struts(index,2),:));
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
        function potentialStruts = removeNormalConStruts(p0,p1,potentialStruts,struts,verts,tol)
            % remove struts that are already joined properly
            numPotential = length(potentialStruts);
            test = false(numPotential,1);
            for inc = 1:numPotential
                index = potentialStruts(inc);
                q0 = verts(struts(index,1),:);
                q1 = verts(struts(index,2),:);
                diff1 = all(abs(p0-q0)<tol);
                diff2 = all(abs(p1-q0)<tol);
                diff3 = all(abs(p0-q1)<tol);
                diff4 = all(abs(p1-q1)<tol);
                test(inc) = any([diff1,diff2,diff3,diff4]);
            end
            potentialStruts(test) = [];
        end
        function potentialStruts = removeParralelStruts(p0,p1,potentialStruts,struts,verts)
            v = (p0-p1)/norm(p0-p1);
            numPotential = length(potentialStruts);
            test = false(numPotential,1);
            for inc = 1:numPotential
                index = potentialStruts(inc);
                q0 = verts(struts(index,1),:);
                q1 = verts(struts(index,2),:);
                u = (q0-q1)/norm(q0-q1);
                test(inc) = all(abs(v-u)<1e-6);
            end
            potentialStruts(test) = [];
        end
        function potentialStruts =  findSplitStrutBB(splitStruts,potentialStruts,currentBB,currentP1,currentP2,strutsOut,vertsOut,tol)
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
                boundingBox = PLG.getBound(vertsOut,strutsOut,desired);
                isIntersect = PLG.findBbIntersect(currentBB,boundingBox);
                potentialNewStruts = desired(isIntersect);
                potentialNewStruts = PLG.removeNormalConStruts(currentP1,currentP2,potentialNewStruts,strutsOut,vertsOut,tol);
                % dont need to remove parallel as these are already removed in whole strut check
                if isempty(potentialNewStruts)
                    % nothing
                else
                    potentialStruts = [potentialStruts;...
                        [potentialNewStruts',ones(length(potentialNewStruts),1)*2]];
                end
                
            end
        end
        function [vertsTmpOut,strutsTmpOut]= findIntersect(p0,p1,potentialStruts,struts,verts,strutsOut,vertsOut,tol)
            % for struts that have a coincident bounding box determine if they
            % intersect
            numPotentialStruts = size(potentialStruts,1);
            coordinates = zeros(numPotentialStruts,3);
            order = zeros(numPotentialStruts,1);
            for inc = 1:numPotentialStruts
                strutType = potentialStruts(inc,2);
                index = potentialStruts(inc,1);
                if strutType==1
                    q0 = verts(struts(index,1),:);
                    q1 = verts(struts(index,2),:);
                else
                    q0 = vertsOut(strutsOut(index,1),:);
                    q1 = vertsOut(strutsOut(index,2),:);
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
            numNewstruts = length(order)+1;
            vertsTmpOut = [p1;coordinates;p0];
            strutsTmpOut = [1:numNewstruts;2:(numNewstruts+1)]';
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
        function [distance, fromP1, coord1, coord2] = distBetween2Segment(p0, p1, q0, q1)
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
