classdef strutOrientationPLG < PLG
    properties
        strutOrientationOffset
    end
    methods
        % refer to the standard PLG for any methods that are not here.
        function obj = load(obj,file)
            % load Loads data for use in PLG as a specialization this
            % method only loads xlsx files.
            %    Returns a modified PLG object.
            parts = strsplit(file,'.');
            extension = parts{end};
            switch extension
                case obj.loadExtensions{5,1}
                    % excel - assumes lattice model
                    warning('PLG:load','A excel will be loaded assuming it is the same format as a .lattice file.');
                    data = xlsread(file);
                otherwise
                    error('PLG:load','Not a suitable load format only a .xlsx file is allowed in strutOrientatioPLG.');
            end
            numNodes=data(1,1);
            numLinks=data(2,1);
            % data
            obj.vertices = data(3:numNodes+2,1:3);
            obj.struts    = data(numNodes+3:numNodes+numLinks+2,1:2);
            obj.strutDiameter = data(numNodes+3:numNodes+numLinks+2,3);
            obj.strutOrientationOffset = data(numNodes+3:numNodes+numLinks+2,4);

            % sphere diameter must be supplied
            obj.sphereDiameter = data(3:numNodes+2,4);
            if any(isnan(obj.sphereDiameter))
                error('PLG:load','With the strutOrientatioPLG extension the sphere diameter must be provided in the input xlsx file.');
            end
        end
    end
    methods (Access=protected)
        function writeSingleStrut(obj,fid,inc)
            % based on a strut increment write a single strut to file 
            point1 = obj.vertices(obj.struts(inc,1),:);
            point2 = obj.vertices(obj.struts(inc,2),:);
            offset_start = obj.strutOrientationOffset(inc);
            
            strut_vector = point2-point1;
            normalized_strut_vector = strut_vector/norm(strut_vector);
            if normalized_strut_vector(3)==1 || normalized_strut_vector(3)==-1 % vector is aligned with z axis
                reference_vector = [1,0,0];
            else
                reference_vector = cross([0,0,1],normalized_strut_vector);
                reference_vector = reference_vector/norm(reference_vector);
            end
            
            % generate all the points on the strut
            offset = obj.strutDiameter(inc)*reference_vector/2; % a offset for the radius
            vert1end = zeros(obj.resolution,3);% end 1
            vert2end = zeros(obj.resolution,3);% end 2
            for j=1:obj.resolution
                Qrot1 = PLG.qGetRotQuaternion((j-1)*2*pi/obj.resolution+offset_start, normalized_strut_vector);
                absolutePointRotation = PLG.qRotatePoint(offset', Qrot1)';
                
                vert1end(j,:)=absolutePointRotation+point1;
                vert2end(j,:)=absolutePointRotation+point2;
            end
            
            % write the points to file in the right order to create a stl file
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
    end

end