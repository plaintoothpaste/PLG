classdef manufacturabilityColour < PLG
    % enables the custom colouring of a lattice file based on manufaturability
    % constraints
    properties
        incline
        strutLength
        colour
        processMap
    end
    methods
        function obj = manufacturabilityColour(file)
            obj = obj@PLG(file);
            
            % set resolution
            obj = set(obj,'resolution',8);
            obj = set(obj,'sphereResolution',8);
        end
        function obj = runManufacturability(obj,file)
            % convert angle incline and diameter to colours
            % file is an input csv document that specifies manufacturing
            % capability
            % anything outside of the process map is red
            % 'r' - red bad
            % 'y' - borderline
            % 'c' - likely good
            % 'g' - good
            
            obj = readProcessMap(obj,file); % load the process map
            obj = calcInclineAndLength(obj); % calculate length and incline
            
            % get a list of colours
             obj = interpAllColour(obj);
             obj.colour = manufacturabilityColour.letter2numColour(obj.colour); % change letter colours to array colours
        end
        
    end
    methods (Access = protected) % sub functions
        function obj = calcInclineAndLength(obj)
            % calculate length and incline for each face
            v1 = obj.vertices(obj.struts(:,1),:);
            v2 = obj.vertices(obj.struts(:,2),:);
            
            % angle
            angleVector = [0,0,1];
            numFaces = size(obj.struts,1);
            v = v1-v2;
            vCell = mat2cell(v,ones(numFaces,1),3);
            costheta = cellfun(@(v) dot(angleVector,v)./(norm(angleVector).*norm(v)),vCell);
            theta = acos(costheta)*180/pi;
            t = theta>90;
            theta(t) = theta(t)-90;
            theta(~t) = 90-theta(~t);
            obj.incline = theta;
            
            % length
            diffSquare = (v2-v1).^2;
            obj.strutLength = sqrt(sum(diffSquare,2));
        end
        function obj = readProcessMap(obj,file)
            % read in a process map csv and create a lookup table that
            % returns a number based on incline diameter and strut length
            fid = fopen(file,'r');
            t = textscan(fid,'%f %f %f %*[^\n]',1,'Delimiter',',');
            numTables = t{1};
            rows = t{2}+1;
            cols = t{3}+1;
            tot = rows*numTables;
            t = textscan(fid,'%s',tot+1,'Delimiter','\n');
            fclose(fid);
            
            obj.processMap = manufacturabilityColour.readTables(t{1}(2:end),numTables,rows,cols);
        end
        function [inclineI,diaI,spanI] = interpIndex(obj,inc)
            % based on the manufacturing options return the correct colour index
            cInc = obj.incline(inc);
            cDia = obj.strutDiameter(inc);
            cSpan = obj.strutLength(inc);
            
            [~,inclineI] = min(abs(obj.processMap.alpha-cInc));
            [~,diaI] = min(abs(obj.processMap.dia-cDia));
            [~,spanI] = min(abs(obj.processMap.span-cSpan));
        end
        function obj = interpAllColour(obj)
            % interp all the colours
            numStruts = length(obj.incline);
            obj.colour = cell(numStruts,1);
            for inc = 1:numStruts
                [inclineI,diaI,spanI] = interpIndex(obj,inc);
                obj.colour{inc} = obj.processMap.colour{spanI,inclineI,diaI};
            end
        end
    end
    methods (Static)
        function d = readTables(data,numTables,rows,cols)
            % splits the raw data from the text file into usefull
            % information
            d.dia = zeros(numTables,1);
            d.alpha = zeros(1,cols-1);
            d.span = zeros(rows-1,1);
            d.colour = cell(rows-1,cols-1,numTables);
            formatSpec = ['%f ',repmat('%s ',1,cols-1)];
            formatSpecHeader = ['%f ',repmat('%f ',1,cols-1)];
            for inc = 1:numTables
                cData = data(1:rows); 
                data(1:rows)=[];
                headerRow = textscan(cData{1},formatSpecHeader,'Delimiter',',');
                d.dia(inc) = headerRow{1};
                normalData = cellfun(@(str) textscan(str,formatSpec,'Delimiter',','),cData(2:end),'UniformOutput',0);
                % apply colours
                for incR = 1:rows-1
                    nd = normalData{incR};
                    for incC = 1:cols-1
                        d.colour{incR,incC,inc} = nd{incC+1}{1};
                    end
                end
                
                if inc==1
                    d.alpha = cell2mat(headerRow(2:end));
                    d.span = cellfun(@(x) x{1},normalData);
                else
                    test1 = any(d.alpha~=cell2mat(headerRow(2:end)));
                    test2 = any(d.span~=cellfun(@(x) x{1},normalData));
                    if test1 || test2
                        error('manufacturabilityColour:inconsistentTables','Input process map tables are inconsistent');
                    end
                end
            end
        end
        function colourNum = letter2numColour(colourLet)
            % convert a list of colours in a cell, as letters, to an array of colour values
            if isempty(colourLet)
                colourNum = [];
                return;
            end
            cLetter = colourLet{1};
            if strcmp(cLetter,'g')
                colourNumC = [0,1,0];
            elseif strcmp(cLetter,'c')
                colourNumC = [0,1,1];
            elseif strcmp(cLetter,'y')
                colourNumC = [1,1,0];
            elseif strcmp(cLetter,'r')
                colourNumC = [1,0,0];
            else
                error('letter %s not supported',cLetter);
            end
            
            colourNumF = manufacturabilityColour.letter2numColour(colourLet(2:end));
            colourNum = [colourNumC;colourNumF];
        end
    end
    methods (Access = ?TestManufacturabilityColour)
        % interfaces for testing only
        function actRes = readProcMap(obj,file)
            obj = readProcessMap(obj,file);
            actRes = obj.processMap;
        end
        function actRes = testLengthIncline(obj,file)
            obj = readProcessMap(obj,file); % load the process map
            obj = calcInclineAndLength(obj); % calculate length and incline
            actRes.dia = obj.strutDiameter;
            actRes.len = obj.strutLength;
        end
        function actRes = testInterpIndex(obj,file,incer)
            obj = readProcessMap(obj,file); % load the process map
            obj = calcInclineAndLength(obj); % calculate length and incline
            
            numInc = length(incer);
            inclineI = zeros(numInc,1);
            diaI = zeros(numInc,1);
            spanI = zeros(numInc,1);
            for inc = 1:numInc
                [inclineI(inc),diaI(inc),spanI(inc)] = interpIndex(obj,incer(inc));
            end
            
            actRes.inclineI = inclineI;
            actRes.diaI = diaI;
            actRes.spanI = spanI;
        end
    end
    
end