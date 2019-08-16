classdef manufacturabilityColour < PLG
    % enables the custom colouring of a file based on manufaturability
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
            obj = calcInclineAndLength(obj);
            obj.colour = zeros(size(obj.incline));
            
            test = obj.incline<30;
            obj.colour(test) = 63488; % red
            
            test = obj.incline>=30 & obj.incline<60;
            obj.colour(test) = 62; % blue
            
            test = obj.incline>=60;
            obj.colour(test) = 1984; % green
        end
    end
    methods (Access = protected) % sub functions
        function obj = calcInclineAndLength(obj)
            angleVector = [0,0,1];
            numFaces = size(obj.struts,1);
            v = obj.vertices(obj.struts(:,1),:)-obj.vertices(obj.struts(:,2),:);
            vCell = mat2cell(v,ones(numFaces,1),3);
            costheta = cellfun(@(v) dot(angleVector,v)./(norm(angleVector).*norm(v)),vCell);
            theta = acos(costheta)*180/pi;
            t = theta>90;
            theta(t) = theta(t)-90;
            theta(~t) = 90-theta(~t);
            obj.incline = theta;
            
            obj.strutLength = '';
        end
        function obj = readProcessMap(obj,file)
            % read in a process map csv and create a lookup table that
            % returns a number based on incline diameter and strut length
            fid = fopen(file,'r');
            t = textscan(fid,'%f %f %f %*[^\n]',1,'Delimiter',',');
            numTables = t{1};
            rows = t{2};
            cols = t{3};
            tot = rows*numTables;
            t = textscan(fid,'%s',tot+1,'Delimiter','\n');
            fclose(fid);
            
            data = manufacturabilityColour.readTables(t{1}(2:end),numTables,rows,cols);
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
                        d.colour{incR,incC,inc} = nd{incC+1};
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
    end
    methods (Access = ?TestManufacturabilityColour)
        % interfaces for testing only
        function expRes = readProcMap(obj,file)
            obj = readProcessMap(obj,file);
            expRes = obj.processMap;
        end
    end
    
end