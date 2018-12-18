classdef colourMod < PLG
    % enables the custom colouring of a stl file based on incline
    properties
        incline
        colour
    end
    methods
        function obj = colourMod(file)
            obj = obj@PLG(file);
            obj = calcIncline(obj);
            obj = setColours(obj);
            
            % set resolution
            obj = set(obj,'resolution',8);
            obj = set(obj,'sphereResolution',8);
            
            fileOut = [file,'.stl'];
            saveStl(obj,fileOut);
        end
    end
    methods (Access=protected) % stl format
        function obj = calcIncline(obj)
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
        end
        function obj = setColours(obj)
            % convert angle to colour
            obj.colour = zeros(size(obj.incline));
            
            test = obj.incline<30;
            obj.colour(test) = 63488; % red
            
            test = obj.incline>=30 & obj.incline<60;
            obj.colour(test) = 62; % blue
            
            test = obj.incline>=60;
            obj.colour(test) = 1984; % green
        end
        function writeSingleFace(obj,fid,v1,v2,v3,inc)
            normal=cross(v2-v1,v3-v1);
            colour = obj.colour(inc);
            fwrite(fid,normal,'float32');           % write normal vector floating point numbers
            fwrite(fid,v1,'float32');   % first vertex (x,y,z)
            fwrite(fid,v2,'float32');   % second vertex
            fwrite(fid,v3,'float32');   % Third vertex
            fwrite(fid,colour,'uint16','l');
        end
    end
end