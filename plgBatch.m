classdef plgBatch < matlab.unittest.TestCase
    % this class enables the batch creation of a large number of lattice
    % structures automatically using either sequential or exhaustive
    % combinations.
    %
    % to run a specific function use:
    % results = run(plgBatch,'functionNameHere');
    %
    % examples include:
    % results = run(plgBatch,'runAllCombinations');
    % results = run(plgBatch,'squareUnitCell');
    % results = run(plgBatch,'squareLattice');
    properties 
        % user set properties
        projectName = 'doeLattice';
        outputFolder = 'C:\scratch\CODE\tmp\batchOut'; %'./batchOut'
        useSpheres = true;
        supportDia = 0.25;
        useBaseFlat = true;
    end
    properties (TestParameter)
        % every combination of properties will be produced must be a cell
        % input
        saveOut = {'custom'};
        unitCell = {{'verticalFaceRods','zRods'}};
        resolution = {15};
        
        unitSizeX = {5}
        unitSizeY = {4};
        unitSizeZ = {[100,25,4,1]};
        
        repsX = {1,2,5,10};
        repsY = {4};
        repsZ ={100,25,4,1};
        
        strut_dia = num2cell(0.5);
        ball_dia = num2cell(0.5);
        
    end
    methods (Test, ParameterCombination='exhaustive')
        function runAllCombinations(obj,saveOut,unitCell,resolution,unitSizeX,unitSizeY,unitSizeZ,repsX,repsY,repsZ,strut_dia,ball_dia)
            % generate the plg object with the desired inputs then save
            % with the desired file formats to the output folder
            
            plgObj = PLG();
            plgObj = set(plgObj,'resolution',resolution);
            plgObj = set(plgObj,'strutDiameter',strut_dia);
            plgObj = set(plgObj,'sphereAddition',obj.useSpheres);
            plgObj = set(plgObj,'sphereDiameter',ball_dia);
            plgObj = set(plgObj,'sphereResolution',resolution);
            plgObj = set(plgObj,'unitSize',[unitSizeX,unitSizeY,unitSizeZ]);
            plgObj = set(plgObj,'origin',[0,0,0]);
            plgObj = set(plgObj,'replications',[repsX,repsY,repsZ]);
            plgObj = set(plgObj,'baseFlat',obj.useBaseFlat);
            plgObj = defineUnit(plgObj,unitCell);
            plgObj = cellReplication(plgObj);
            plgObj = cleanLattice(plgObj);
            if obj.useSpheres
                fileName = sprintf('unit_%s d_%04.1f usx_%04.1f usy_%04.1f usz_%04.1f rx_%04.1f ry_%04.1f rz_%04.1f_bd_%04.1f',plgObj.unitName,strut_dia,unitSizeX,unitSizeY,unitSizeZ,repsX,repsY,repsZ,ball_dia);
            else
                fileName = sprintf('unit_%s d_%04.1f usx_%04.1f usy_%04.1f usz_%04.1f rx_%04.1f ry_%04.1f rz_%04.1f',plgObj.unitName,strut_dia,unitSizeX,unitSizeY,unitSizeZ,repsX,repsY,repsZ);
            end
            
            
            fileLocation = sprintf('%s%s%s',obj.outputFolder,filesep,fileName);
            switch saveOut
                case 'stl'
                    saveStl(plgObj,[fileLocation,'.stl']);
                case 'custom'
                    saveCustom(plgObj,[fileLocation,'.custom']);
                case '3mf'
                    save3mf(plgObj,[fileLocation,'.3mf']);
                case 'all'
                    saveStl(plgObj,[fileLocation,'.stl']);
                    saveCustom(plgObj,[fileLocation,'.custom']);
                    save3mf(plgObj,[fileLocation,'.3mf']);
                otherwise
                    error('saveOut type:%s not supported',saveOut);
            end
        end
        function squareUnitCell(obj,saveOut,unitCell,resolution,unitSizeX,repsX,repsY,repsZ,strut_dia,ball_dia)
            % as a square unit cell is quite common this option is the same
            % as the above but with only a unitSizeX used for all
            % dimensions
            plgObj = PLG();
            plgObj = set(plgObj,'resolution',resolution);
            plgObj = set(plgObj,'strutDiameter',strut_dia);
            plgObj = set(plgObj,'sphereAddition',obj.useSpheres);
            plgObj = set(plgObj,'sphereDiameter',ball_dia);
            plgObj = set(plgObj,'sphereResolution',resolution);
            plgObj = set(plgObj,'unitSize',[unitSizeX,unitSizeX,unitSizeX]);
            plgObj = set(plgObj,'origin',[0,0,0]);
            plgObj = set(plgObj,'replications',[repsX,repsY,repsZ]);
            plgObj = set(plgObj,'baseFlat',obj.useBaseFlat);
            plgObj = defineUnit(plgObj,unitCell);
            plgObj = cellReplication(plgObj);
            plgObj = cleanLattice(plgObj);
            
            fileName = sprintf('unit_%s d_%04.1f us_%04.1f rx_%04.1f ry_%04.1f rz_%04.1f',plgObj.unitName,strut_dia,unitSizeX,repsX,repsY,repsZ);
            fileLocation = sprintf('%s%s%s',obj.outputFolder,filesep,fileName);
            switch saveOut
                case 'stl'
                    saveStl(plgObj,[fileLocation,'.stl']);
                case 'custom'
                    saveCustom(plgObj,[fileLocation,'.custom']);
                case '3mf'
                    save3mf(plgObj,[fileLocation,'.3mf']);
                case 'all'
                    saveStl(plgObj,[fileLocation,'.stl']);
                    saveCustom(plgObj,[fileLocation,'.custom']);
                    save3mf(plgObj,[fileLocation,'.3mf']);
                otherwise
                    error('saveOut type:%s not supported',saveOut);
            end
        end
        function squareLattice(obj,saveOut,unitCell,resolution,unitSizeX,repsX,strut_dia)
            % ball diameter will just match strut
            plgObj = PLG();
            plgObj = set(plgObj,'resolution',resolution);
            plgObj = set(plgObj,'strutDiameter',strut_dia);
            plgObj = set(plgObj,'sphereAddition',obj.useSpheres);
            plgObj = set(plgObj,'sphereDiameter',strut_dia);
            plgObj = set(plgObj,'sphereResolution',resolution);
            plgObj = set(plgObj,'unitSize',[unitSizeX,unitSizeX,unitSizeX]);
            plgObj = set(plgObj,'origin',[0,0,0]);
            plgObj = set(plgObj,'replications',[repsX,repsX,repsX]);
            plgObj = set(plgObj,'baseFlat',obj.useBaseFlat);
            plgObj = defineUnit(plgObj,unitCell);
            plgObj = cellReplication(plgObj);
            plgObj = cleanLattice(plgObj);
            
            fileName = sprintf('unit_%s d_%04.1f us_%04.1f rep_%04.1f',plgObj.unitName,strut_dia,unitSizeX,repsX);
            fileLocation = sprintf('%s%s%s',obj.outputFolder,filesep,fileName);
            switch saveOut
                case 'stl'
                    saveStl(plgObj,[fileLocation,'.stl']);
                case 'custom'
                    saveCustom(plgObj,[fileLocation,'.custom']);
                case '3mf'
                    save3mf(plgObj,[fileLocation,'.3mf']);
                case 'all'
                    saveStl(plgObj,[fileLocation,'.stl']);
                    saveCustom(plgObj,[fileLocation,'.custom']);
                    save3mf(plgObj,[fileLocation,'.3mf']);
                otherwise
                    error('saveOut type:%s not supported',saveOut);
            end
        end
        function simulationCompare(obj,saveOut,unitCell,resolution,unitSizeX,repsX,repsZ,strut_dia)
            % modified version of square lattice that adds support
            plgObj = PLG();
            plgObj = set(plgObj,'resolution',resolution);
            plgObj = set(plgObj,'strutDiameter',strut_dia);
            plgObj = set(plgObj,'sphereAddition',obj.useSpheres);
            plgObj = set(plgObj,'sphereDiameter',strut_dia);
            plgObj = set(plgObj,'sphereResolution',resolution);
            plgObj = set(plgObj,'unitSize',[unitSizeX,unitSizeX,unitSizeX]);
            plgObj = set(plgObj,'origin',[0,0,0]);
            plgObj = set(plgObj,'replications',[repsX,repsX,repsZ]);
            plgObj = set(plgObj,'baseFlat',obj.useBaseFlat);
            plgObj = defineUnit(plgObj,unitCell);
            plgObj = cellReplication(plgObj);
            plgObj = cleanLattice(plgObj);
            
            fileName = sprintf('unit_%s d_%04.1f us_%04.1f xy_%04.1f z_%04.1f ',plgObj.unitName,strut_dia,unitSizeX,repsX,repsZ);
            fileLocation = sprintf('%s%s%s',obj.outputFolder,filesep,fileName);
            saveCustom(plgObj,[fileLocation,'.custom']);
            
            plgObj = addSupport([fileLocation,'.custom'],obj.supportDia,0,10,0.1);
            plgObj = padSupport(plgObj,0.9+strut_dia/2,obj.supportDia,0);
            plgObj = set(plgObj,'sphereResolution',resolution);
            plgObj = set(plgObj,'resolution',resolution);
            
            switch saveOut
                case 'stl'
                    saveStl(plgObj,[fileLocation,'.stl']);
                case 'custom'
                    saveCustom(plgObj,[fileLocation,'.custom']);
                case '3mf'
                    save3mf(plgObj,[fileLocation,'.3mf']);
                case 'all'
                    saveStl(plgObj,[fileLocation,'.stl']);
                    saveCustom(plgObj,[fileLocation,'.custom']);
                    save3mf(plgObj,[fileLocation,'.3mf']);
                otherwise
                    error('saveOut type:%s not supported',saveOut);
            end
        end
    end
    methods (Test, ParameterCombination='sequential')
        % functions in this method block will only run if the inputs are
        % all the same length. and each input will be paired
        
    end
end