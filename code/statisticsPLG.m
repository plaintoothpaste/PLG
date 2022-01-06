classdef statisticsPLG < manufacturablePLG
    properties
        maxwell_number
        aspect
        status
    end
    methods
        function obj = statisticsPLG(lattice_file,process_map_file)
            % first run manufacturability
            obj = obj@manufacturablePLG(lattice_file);
            obj = readProcessMap(obj,process_map_file);
            obj = cleanLattice(obj);
            obj = calcInclineAndLength(obj);
            obj = interpAllColour(obj); % calculates an array of letters for each strut under obj.colour
            
            % next run extra calculations
            obj = calcMaxwell(obj);
            obj.aspect = obj.strutLength./obj.strutDiameter;
            obj = colour2status(obj);
        end
        function save(obj,file_name)
            % save exports the properties to a xlsx file
            %    This method overloads thhe default save behaviour of PLG
            
            % Maxwell to output format
            if obj.maxwell_number>0
                LatBehav='Over-stiff';
            elseif obj.maxwell_number==0
                LatBehav='Just-stiff';
            elseif obj.maxwell_number<0
                LatBehav='Under-stiff';
            end
            
            % counters of strut types
            n_failed = sum(cellfun(@(s) strcmp(s,'Failed zone'),obj.status,'UniformOutput',true));
            n_comp = sum(cellfun(@(s) strcmp(s,'Compromised zone'),obj.status,'UniformOutput',true));
            n_robust = sum(cellfun(@(s) strcmp(s,'Robust zone'),obj.status,'UniformOutput',true));
            n_over_robust = sum(cellfun(@(s) strcmp(s,'Over robust zone'),obj.status,'UniformOutput',true));

            summary_sheet = {
                'Maxwell Number = ',obj.maxwell_number;
                'Maxwell Type = ',LatBehav;
                'Number of Failed struts = ',n_failed;
                'Number of Compromised struts = ',n_comp;
                'Number of Robust struts = ',n_robust;
                'Number of Over robust struts = ',n_over_robust;
                };

            %% create a page for the properties of every strut
            header={'Strut ID','Node 1','Node 2','Strut Diameter','Strut Length','Angle to Platen','Strut Quality Manufacture Zone','Aspect Ratio'};
            n = size(obj.struts,1);
            data = cell(n,8);
            data(:,1) = num2cell([1:n]');
            data(:,2) = num2cell(obj.struts(:,1));
            data(:,3) = num2cell(obj.struts(:,2));
            data(:,4) = num2cell(obj.strutDiameter);
            data(:,5) = num2cell(obj.strutLength);
            data(:,6) = num2cell(obj.incline);
            data(:,7) = obj.status;
            data(:,8) = num2cell(obj.aspect);
            
            %% write all the data
            xlswrite(file_name,summary_sheet,'Sheet1'); % default sheet
            xlswrite(file_name,header,'All_struts');
            xlswrite(file_name,data,'All_struts','A2');
        end
    end

    methods (Access = private)
        function obj = calcMaxwell(obj)
            % Checks if nodes lie on a single plane and calculate maxwell
            numStruts = size(obj.struts,1);
            numNodes = size(obj.vertices,1);
            nodes=[obj.vertices(:,1), obj.vertices(:,3), obj.vertices(:,2), obj.sphereDiameter];
            if length(unique(nodes(:,1)))==1 || length(unique(nodes(:,2)))==1 || length(unique(nodes(:,3)))==1
                obj.maxwell_number=numStruts-(2*numNodes)+3;
            else
                % 3D
                obj.maxwell_number=numStruts-(3*numNodes)+ 6;
            end
        end
        function obj = colour2status(obj)
            % convert the colour codes in manufacturable PLG to words for
            % exporting to excel
            %     'r' - Failed zone
            %     'y' - Compromised zone
            %     'c' - Robust zone
            %     'g' - Over robust zone
            obj.status = cellfun(@(c) statisticsPLG.c2w(c),obj.colour,'UniformOutput',false);
        end
    end
    methods (Static)
        function word = c2w(colour)
            switch colour
                case "r"
                    word = "Failed zone";
                case "y"
                    word = "Compromised zone";
                case "c"
                    word = "Robust zone";
                case "g"
                    word = "Over robust zone";
                otherwise
                    error('letter %s not supported',colour);
            end
        end
    end
end



