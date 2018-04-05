classdef PLG
    % The PLG program designed for the generation of regular repeating stuctures, normally a lattice
    % but not limited to lattices.
    properties (SetAccess=protected)
        unitName;
        unitType; % beam - transform storage of single beam can be saved to all formats 
                  % facet - all facets can save out as stl or 3mf but 3mf not recommended
                  % unit - one or more unit cells defined in facet format and
                  % then transformed to the desired location 3mf or stl
                  % custom - data from a custom beam input file can be
                  % converted to beam for save out
        
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
        
        loadExtensions = {'xml','custom unit cell file defined as a xml';...
                          'stl', 'standard unit cell file';...
                          'custom','standard beam output of PLG';...
                          'csv','manualy defined beam model';...
                          'xlsx','manualy defined beam model'};
        saveExtensions = {'stl', 'binary facet representation (compatibility)';...
                          '3mf', 'tesselated facet file (recommended)';...
                          'custom', 'beam output method (simulation)';...
                          'inp', 'Abaqus input file (TODO)'};
                      
        tolerance; % defined as 1/100 of the shortest length present
        dx; % length of a strut
    end
    methods
        function obj = PLG(varargin)
            % creates the PLG object
            switch numel(varargin)
                case 0
                    % generate a new lattice
                    disp('generating a lattice from scratch');
                    obj.sphereAddition = false;
                    obj.unitType = 'beam';
                case 1
                    % import a custom lattice file containing beam and node
                    % definitions see load function for more information
                    disp('Loading an existing file');
                    obj = load(obj,varargin{1});
                    obj.unitType = 'custom';
                    obj.sphereAddition = true;
                    obj.unitName = varargin{1};
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
        function obj = defineUnit(obj,unitNames,type)
            % define a unit cell and wheter output should be beam(default)
            % unit or facet
            if ~exist('type','var')
                type = 'beam';
            end
            
            addpath('unitCell');
            unitObj = unitCell(unitNames,obj);
            obj = scale(obj);
            switch type
                case 'beam'
                    [vertices,connections,transform,name,type] = beamOut(unitObj);
                case 'unit'
                    [vertices,connections,transform,name,type] = unitOut(unitObj);
                case 'facet'
                    [vertices,connections,transform,name,type] = facetOut(unitObj);
                otherwise
                    error('not a suitable type');
            end
            rmpath('unitCell');
            
            % place onto object
            obj.unitType = type;
            obj.transform = transform;
            obj.vertices = vertices;
            obj.struts = connections;
            obj.unitName = name;
        end
        function obj = cellReplication(obj)
            % if unitType is a beam it will replicate transformation if it is a facet type it will
            % replicate struts and verts
            % replaces obj.replications with a list of x,y,z coordinates for the unit cell
            xPlacement = [0:obj.unitSize(1):obj.unitSize(1)*(obj.replications(1)-1)]+obj.origin(1);
            yPlacement = [0:obj.unitSize(2):obj.unitSize(2)*(obj.replications(2)-1)]+obj.origin(2);
            zPlacement = [0:obj.unitSize(3):obj.unitSize(3)*(obj.replications(3)-1)]+obj.origin(3);
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
                case 'custom'
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
                    error('cell replication can not be performed with this method or a unit cell is not yet defined');
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
                case 'custom'
                    % tolerance - get the shortest strut and divide by 1000
                    verts1 = obj.vertices(obj.struts(:,1),:);
                    verts2 = obj.vertices(obj.struts(:,2),:);
                    lengthVerts = sum(sqrt((verts1-verts2).^2),2);
                    obj.tolerance = min(lengthVerts)/15;
                    % duplicate vertices
                    [obj.vertices,i,indexn]=uniquetol(obj.vertices,obj.tolerance,'ByRows',1,'DataScale',1);
                    if ~isempty(obj.sphereDiameter)
                        if numel(obj.sphereDiameter)==1
                            obj.sphereDiameter = ones(size(verts1,1),1)*obj.sphereDiameter;
                        end
                        obj.sphereDiameter = obj.sphereDiameter(i);
                    end
                    obj.struts = indexn(obj.struts);
                    
                    % duplicate struts
                    tmp = obj.struts;
                    test = tmp(:,1)>tmp(:,2);
                    obj.struts(test,1) = tmp(test,2);
                    obj.struts(test,2) = tmp(test,1);
                    [obj.struts,i] = unique(obj.struts,'rows');
                    if numel(obj.strutDiameter)==1
                        obj.strutDiameter = ones(size(tmp,1),1)*obj.strutDiameter;
                    end
                    obj.strutDiameter = obj.strutDiameter(i);
                    
                    % duplicate struts zero length
                    test = obj.struts(:,1)==obj.struts(:,2);
                    obj.struts(test,:)=[];
                    obj.strutDiameter(test)=[];
                otherwise
                    error('cleaning can not be performed on this unitType');
            end
        end
        
        function obj = beam2custom(obj)
            % changes a beam file to the custom import format
            % this enables stl and custom save
            numTransforms = size(obj.transform,1);
            newVerts = zeros(numTransforms*2,3);
            newFacets = [(1:2:(numTransforms*2-1))',(2:2:numTransforms*2)'];
            locations = [obj.vertices,[1;1]];
            for inc = 1:numTransforms
                squareTransform = flat2squareTransform(obj,inc);
                newLocation = locations*squareTransform;
                newVerts((2*inc-1):2*inc,:) = newLocation(:,1:3);
            end
            % replace vertices and facets change unit type and clear transforms
            obj.unitType = 'custom';
            obj.vertices = newVerts;
            obj.struts = newFacets;
            obj.transform = [];
        end
        function obj = custom2beam(obj)
            % converts a custom loaded file to a beam format which stores data as a series of
            % transforms of the generic beam file
            
            % as this is already completed in the unit cell function simply use it
            addpath('unitCell');
            unitObj = unitCell({'bcc'},obj);
            unitObj.vertices = obj.vertices;
            unitObj.connections = obj.struts;
            
            [vertices,connections,transform,~,~] = beamOut(unitObj);
            rmpath('unitCell');
            
            obj.unitType = 'beam';
            obj.transform = transform;
            obj.vertices = vertices;
            obj.struts = connections;
        end
        function obj = beam2facet(obj)
            % not recommended but kept for compatibility represents all the facets in the 3d
            % geometry all the time.
            error('TODO');
        end
        function obj = custom2facet(obj)
            error('TODO');
        end
        function obj = unit2beam(obj)
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
                case 'custom'
                    f = figure;
                    f.Units	= 'Normalized';
                    f.Position = [0,0,1,1];
                    f.Name = 'STL plotter';
                    a = axes;
                    a.View = [45,45];
                    axis vis3d
                    a.NextPlot='add';
                    numStruts = length(obj.strutDiameter);
                    if ~exist('colours','var')
                        colours = repmat([0.3,0.3,0.3,0.5],numStruts,1);
                    end
                    for inc = 1:numStruts
                        points = [obj.vertices(obj.struts(inc,1),:);obj.vertices(obj.struts(inc,2),:)];
                        p = plot3(points(:,1),points(:,2),points(:,3),'Color',colours(inc,:));
                        p.MarkerFaceColor = [0.9,0.5,0]; p.MarkerEdgeColor = 'none';
                    end
                otherwise
                    error('Plotting not yet supported for this unitType');
            end
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
    end
    methods % lattice manipulations
        function obj = translate(obj,x,y,z)
            %translates a lattice in space
            switch obj.unitType
                case 'beam'
                    obj = beamTranslate(obj,x,y,z);
                case 'custom'
                    obj = standardTranslate(obj,x,y,z);
                case 'facet'
                    obj = standardTranslate(obj,x,y,z);
            end
        end
        function obj = rotate(obj,wx,wy,wz)
            % rotations are in degrees about the main axes
            % to radians
            wx = wx*pi/180;
            wy = wy*pi/180;
            wz = wz*pi/180;
            switch obj.unitType
                case 'beam'
                    obj = beamRotate(obj,wx,wy,wz);
                case 'custom'
                    
                case 'facet'
                    
            end
        end
        function obj = plus(obj,obj1)
            % add obj1 to obj2
            if ~strcmp(obj.unitType,obj1.unitType)
                error('Both PLG objects must be the samae type to add together');
            end
            switch obj.unitType
                case 'beam'
                    % put the transforms together
                    obj.transform = [obj.transform;obj1.transform];
                case 'unit'
                    error('TODO');
                otherwise
                    newStruts = obj1.struts+max(obj.struts(:));
                    obj.struts = [obj.struts;newStruts];
                    obj.strutDiameter = [obj.strutDiameter;obj1.strutDiamter];
                    obj.sphereDiameter = [obj.sphereDiameter;obj1.sphereDiameter];
                    obj.vertices = [obj.vertices;obj1.vertices];
            end
            
            % remove any matching struts/transforms
            obj = cleanLattice(obj);
        end
    end
    methods % stats and advanced
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
        function obj = load(obj,file)
            % load a custom beam input file to generate a lattice structure
            parts = strsplit(file,'.');
            extension = parts{end};
            switch extension
                case obj.loadExtensions{1,1}
                    % xml - use unitCell class
                    %TODO
                    error('use unit generate to load file')
                case obj.loadExtensions{2,1}
                    % stl - use stlHandler class
                    %TODO
                    error('use unitCell/stl2unitCell.m to convert to a xml for loading')
                case obj.loadExtensions{3,1}
                    % custom - assumes custom model type
                    data = csvread(file);
                case obj.loadExtensions{4,1}
                    % csv - assumes custom model type
                    data = csvread(file);
                case obj.loadExtensions{5,1}
                    % exvel - assumes custom model
                    data = xlsread(file);
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
        function obj = standardTranslate(obj,x,y,z)
            obj.vertices(:,1) = obj.vertices(:,1)+x;
            obj.vertices(:,2) = obj.vertices(:,2)+y;
            obj.vertices(:,3) = obj.vertices(:,3)+z;
        end
        function obj = beamTranslate(obj,x,y,z)
            obj.transform(:,10) = obj.transform(:,10)+x;
            obj.transform(:,11) = obj.transform(:,11)+y;
            obj.transform(:,12) = obj.transform(:,12)+z;
        end
        function obj = beamRotate(obj,wx,wy,wz)
            transX = [1,0       ,0      ,0;...
                      0,cos(wx) ,sin(wx),0;...
                      0,-sin(wx),cos(wx),0;...
                      0,0       ,0      ,1];
            transY = [cos(wy),0,-sin(wy),0;...
                      0      ,1,0       ,0;...
                      sin(wy),0,cos(wy) ,0;...
                      0      ,0,0       ,1];
            transZ = [cos(wz) ,sin(wz),0,0;...
                      -sin(wz),cos(wz),0,0;...
                      0       ,0      ,1,0;...
                      0       ,0      ,0,1];
            fullTrans = transX * transY * transZ;
            for inc = 1:size(obj.transform,1)
                squareTransform = flat2squareTransform(obj,inc);
                newTrans = fullTrans*squareTransform;
                obj = squareIntoFlatTransform(obj,newTrans,inc);
            end
        end
        function squareTransform = flat2squareTransform(obj,index)
            % takes the transform at row index and returns the square transform ready for use
            currentTransformFlat = obj.transform(index,:);
            squareTransform = eye(4);
            squareTransform(1,1:3)=currentTransformFlat(1:3);
            squareTransform(2,1:3)=currentTransformFlat(4:6);
            squareTransform(3,1:3)=currentTransformFlat(7:9);
            squareTransform(4,1:3)=currentTransformFlat(10:12);
        end
        function obj = squareIntoFlatTransform(obj,squareTransform,index)
            % takes a square transform and flattens it then replaces the transform at the index row
            obj.transform(index,:) = [squareTransform(1,1:3),...
                squareTransform(2,1:3),...
                squareTransform(3,1:3),...
                squareTransform(4,1:3)];
        end
    end
                
    methods % save out methods
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
                    saveCustom(obj,file)
                case 4
                    saveAbaqus(obj,file)
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
        function saveCustom(obj,fullName)
            % saves data to a csv file with the .custom extension usefull as a beam model input to
            % simulations
            switch obj.unitType
                case 'beam'
                    % convert to custom
                    obj = beam2custom(obj);
                case 'custom'
                    % ready to save
                otherwise
                    error('can not save out with this unit type: %s',obj.unitType);
            end
            obj = cleanLattice(obj);
            
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
            idStrut = 0;
            idStrutComponents = 1;
            threeMfDoc = threeMfUnitObj(obj,threeMfDoc,resourcesNode,idStrut);
            threeMfDoc = threeMfreplicate_unit(obj,threeMfDoc,resourcesNode,idStrut,idStrutComponents);
            % write each object (with all its components 1 time
            buildNode = threeMfDoc.createElement('build');
            threeMfNode.appendChild(buildNode);
            itemNode = threeMfDoc.createElement('item');
            itemNode.setAttribute('objectid',num2str(idStrutComponents));
            buildNode.appendChild(itemNode);
            if obj.sphereAddition
                idBall = 2;
                idBallComponents = 3;
                threeMfDoc = threeMfBallObj(obj,threeMfDoc,resourcesNode,idBall);
                threeMfDoc = threeMfreplicate_Ball(obj,threeMfDoc,resourcesNode,idBall,idBallComponents);
                itemNode = threeMfDoc.createElement('item');
                itemNode.setAttribute('objectid',num2str(idBallComponents));
                buildNode.appendChild(itemNode);
            end
            
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
        function threeMfDoc = threeMfreplicate_unit(obj,threeMfDoc,resourcesNode,idRef,newId)
            % for every transform write it to the file in a object that holds every component of the
            % struts
            
            % setup the object
            objNode = threeMfDoc.createElement('object');
            objNode.setAttribute('id',num2str(newId));
            objNode.setAttribute('type','model');
            resourcesNode.appendChild(objNode);
            
            componentsNode = threeMfDoc.createElement('components');
            objNode.appendChild(componentsNode);
            
            numTransform = size(obj.transform,1);
            for inc = 1:numTransform
                componentNode = threeMfDoc.createElement('component');
                componentNode.setAttribute('objectid',num2str(idRef));
                
                currentTransform = obj.transform(inc,:);
                str = sprintf('%5.5f ',currentTransform);
                str(end) = [];
                componentNode.setAttribute('transform',str);
                componentsNode.appendChild(componentNode);
            end
            
            if ~obj.sphereAddition
                return
            end
        end
        function threeMfDoc = threeMfreplicate_Ball(obj,threeMfDoc,resourcesNode,idRef,newId)
            % write the ball locations
            
            % setup the object
            objNode = threeMfDoc.createElement('object');
            objNode.setAttribute('id',num2str(newId));
            objNode.setAttribute('type','model');
            resourcesNode.appendChild(objNode);
            
            componentsNode = threeMfDoc.createElement('components');
            objNode.appendChild(componentsNode);
            
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
                componentNode = threeMfDoc.createElement('component');
                componentNode.setAttribute('objectid',num2str(idRef));
                
                str = sprintf('1 0 0 0 1 0 0 0 1 %5.5f %5.5f %5.5f', verts(inc,1), verts(inc,2),verts(inc,3));
                componentNode.setAttribute('transform',str);
                componentsNode.appendChild(componentNode);
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
    
end
