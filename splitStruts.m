classdef splitStruts < PLG
    % enables splitting of a bad custom file where beam do interesect but this is not present in the
    % file. splitStruts will identify these and split the beams into 2
    
    properties
    end
    
    methods
        function obj = splitStruts(file)
            % takes a badly defined input lattice and scoures it for any struts that intersect and
            % but do not connect to each other. It then breaks these struts and creates a new clean
            % lattice
            % this functions should only be used if the user knows what they are doing and may
            % destroy unique diameter data
            
            obj = obj@PLG(file);
            if ~strcmp(obj.unitType,'custom')
                error('splitting can only be performed on custom files')
            end
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
    end
    methods (Static)
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

