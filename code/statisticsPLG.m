classdef statisticsPLG < manufacturablePLG & PLG
    properties
        maxwell_number

    end
    methods
        function obj = statisticsPLG(lattice_file,process_map_file)
            % first run manufacturability
            obj = obj@PLG(lattice_file);
            obj = runManufacturability(obj,process_map_file);

            % next calculate properties
            obj = calcProperties(obj);
        end
        function obj = calcProperties(obj)
            % getProperties Perform a statistical analysis.
            % convert plg inputs into the function inputs

            numNodes=size(obj.vertices,1);
            numStruts=size(obj.struts,1); % Read Number of Struts
            Struts=[obj.struts,obj.strutDiameter];
            Nodes=[obj.vertices(:,1) obj.vertices(:,3) obj.vertices(:,2) obj.sphereDiameter];
            
            %% CHECK LATTICE DIMENSIONALITY
            % Checks if nodes lie on a single plane and calculate maxwell
            if length(unique(Nodes(:,1)))==1 || length(unique(Nodes(:,2)))==1 || length(unique(Nodes(:,3)))==1
                obj.maxwell_number=numStruts-(2*numNodes)+3;
            else
                % 3D
                obj.maxwell_number=numStruts-(3*numNodes)+ 6;
            end
            
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
            
        end
        function save(file_name)
            % save exports the properties to a xlsx file
            %    This method overloads thhe default save behaviour of PLG

            % Write for Output Excel Sheet format
            if obj.maxwell_number>0
                LatBehav='Over-stiff';
            elseif obj.maxwell_number==0
                LatBehav='Just-stiff';
            elseif obj.maxwell_number<0
                LatBehav='Under-stiff';
            end
            LatBehavOut=sprintf('%d (%s)',obj.maxwell_number,LatBehav);

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
end



