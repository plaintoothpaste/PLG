classdef addSupport < PLG
    %ADDSUPPORT subclass of PLG for adding support elements.
    %    A custom file can eb loaded and this class will add elements of the users desired diameter from a supported node to the lowest z normal node
    
    properties
        % see PLG
        angleVector = [0,0,1]; % vertical
        supportVertsType % list of vertices that are used for support as well as their type, 1 for a support specific node 2 for a node that is already at zmin
    end
    
    methods
        function obj = addSupport(file,diaStrut,diaSphere,incline,penetrationPercentage)
            % initiate the PLG class
            obj = obj@PLG(file);
            % find nodes that are to be supported
            vertsNeedSupport = smartBaseDectection(obj,incline,penetrationPercentage);
            
            % add support elements
            obj = baseAddition(obj,diaStrut,diaSphere,vertsNeedSupport);
        end
        function obj = padSupport(obj,pad,diaStrut,diaSphere)
            % adds a set distance to all nodes at zmin similiar to real pin supports
            minZ = min(obj.vertices(:,3));
            newMinZ = minZ-pad;
            
            moveVerts = find(obj.supportVertsType==1);
            for inc = moveVerts'
                obj.vertices(moveVerts,3) = newMinZ;
            end
            
            canNotMoveVerts = find(obj.supportVertsType==2);
            numCanNotMoveVerts = length(canNotMoveVerts);
            for inc = 1:numCanNotMoveVerts
                attach = obj.vertices(canNotMoveVerts(inc),:); attach(3) = newMinZ;
                obj.vertices(end+1,:) = attach;
                obj.sphereDiameter(end+1) = diaSphere;
                obj.struts(end+1,:) = [canNotMoveVerts(inc),size(obj.vertices,1)];
                obj.strutDiameter(end+1) = diaStrut;
            end
        end
        function plot(obj)
            % plot with the colours of elements different for support
            colours = repmat([0.8,0,0.2],length(obj.strutDiameter),1);
            supportVerts = find(obj.supportVertsType);
            test = arrayfun(@(x) any(x==supportVerts),obj.struts(:,1)) | arrayfun(@(x) any(x==supportVerts),obj.struts(:,2));
            colours(test,:) = repmat([0,0.6,0],sum(test),1);
            plot@PLG(obj,colours);
            
            % has a bug where the original vertice(s) at zmin will also have all their struts
            % coloured
        end
    end
    methods (Access=protected)
        function vertsNeedSupport = smartBaseDectection(obj,incline,pen)
            minZ = min(obj.vertices(:,3));
            maxZ = max(obj.vertices(:,3));
            plusNpercent = (maxZ-minZ)*pen+minZ;
            inclineLimit = incline*pi/180; % convert setupObj value to radians
            potentialVertsIndex = find(obj.vertices(:,3)<plusNpercent);
            vertsNeedSupport = [];
            for inc = potentialVertsIndex'
                pointCurrent = obj.vertices(inc,:);
                test = obj.struts(:,1)==inc;
                pointsOther = obj.vertices(obj.struts(test,2),:);
                test = obj.struts(:,2)==inc;
                pointsOther = [pointsOther;obj.vertices(obj.struts(test,1),:)];
                
                % calculate vector
                v = [pointsOther(:,1)-pointCurrent(1),pointsOther(:,2)-pointCurrent(2),pointsOther(:,3)-pointCurrent(3)];
                % v = [pointCurrent(1)-pointsOther(:,1),pointCurrent(2)-pointsOther(:,2),pointCurrent(3)-pointsOther(:,3)];
                vCell = mat2cell(v,ones(size(pointsOther,1),1),3);
                costheta = cellfun(@(v) dot(obj.angleVector,v)./(norm(obj.angleVector).*norm(v)),vCell);
                if all(costheta>-inclineLimit) % steeper then incline limit
                    vertsNeedSupport = [vertsNeedSupport,inc];
                end
            end
        end
        function obj = baseAddition(obj,diaStrut,diaSphere,vertsNeedSupport)
            % add elements to the lower nodes that require them and replace supportVertsInd with the
            % list of nodes at zmin that can be extended
            numsupportVertsInd = length(vertsNeedSupport);
            
            % add vertex from each node to z=zmin
            minZ = min(obj.vertices(:,3));
            obj.supportVertsType = zeros(length(obj.sphereDiameter),1);
            for inc = 1:numsupportVertsInd
                currentVert = obj.vertices(vertsNeedSupport(inc),:);
                if abs(currentVert(3)-minZ)<1e-5
                    % close enough do not add a support element
                    obj.vertices(vertsNeedSupport(inc),3)=minZ;
                    obj.supportVertsType(vertsNeedSupport(inc)) = 2;
                else
                    % add a support element and vertex
                    lengthVerts = size(obj.vertices,1)+1;
                    obj.supportVertsType(lengthVerts,1) = 1;
                    
                    attach = currentVert; attach(3) = minZ;
                    obj.vertices(lengthVerts,:) = attach;
                    
                    obj.struts(end+1,:) = [vertsNeedSupport(inc),lengthVerts];
                    
                    obj.sphereDiameter(end+1) = diaSphere;
                    obj.strutDiameter(end+1) = diaStrut;
                end
            end
        end
    end
end