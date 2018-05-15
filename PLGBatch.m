classdef plgBatch < matlab.unittest.TestCase
    properties 
        % user set properties
        projectName = 'doeLattice';
        outputFolder = './batchOut'
        useSpheres = true;
    end
    properties (TestParameter)
        % every combination of properties will be produced must be a cell
        % input
        saveOut = {'stl','custom','3mf'}
        unitCell = {{'verticalFaceRods'},{'verticalFaceRods','zRods'},...
            {'centreCross'},{'centreCross','zRods'}};
        resolution = {10};
        
        unitSizeX = {5,7.5,10};
        unitSizeY = {5,7.5,10};
        unitSizeZ = {5,7.5,10};
        
        repsX = num2cell(1:5);
        repsY = num2cell(1:5);
        repsZ = num2cell(1:5);
        
        strut_dia = num2cell(0.5:0.5:2);
        ball_dia = num2cell(0.5:0.5:2);
        
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
            
            plgObj = defineUnit(plgObj,unitCell);
            plgObj = cellReplication(plgObj);
            plgObj = cleanLattice(plgObj);
            
            fileName = sprintf('unit_%s d_%04.1f usx_%04.1f usy_%04.1f usz_%04.1f rx_%04.1f ry_%04.1f rz_%04.1f',plgObj.unitName,strut_dia,unitSizeX,unitSizeY,unitSizeZ,repsX,repsY,repsZ);
            fileLocation = sprintf('%s%s%s',obj.outputFolder,filesep,fileName);
            switch saveOut
                case 'stl'
                    saveStl(plgObj,[fileLocation,'.stl']);
                case 'custom'
                    saveCustom(plgObj,[fileLocation,'.custom']);
                case '3mf'
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
                otherwise
                    error('saveOut type:%s not supported',saveOut);
            end
        end
        function squareLattice(obj,saveOut,unitCell,resolution,unitSizeX,repsX,strut_dia)
            plgObj = PLG();
            plgObj = set(plgObj,'resolution',resolution);
            plgObj = set(plgObj,'strutDiameter',strut_dia);
            plgObj = set(plgObj,'sphereAddition',obj.useSpheres);
            plgObj = set(plgObj,'sphereDiameter',strut_dia);
            plgObj = set(plgObj,'sphereResolution',resolution);
            plgObj = set(plgObj,'unitSize',[unitSizeX,unitSizeX,unitSizeX]);
            plgObj = set(plgObj,'origin',[0,0,0]);
            plgObj = set(plgObj,'replications',[repsX,repsX,repsX]);
            
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
                otherwise
                    error('saveOut type:%s not supported',saveOut);
            end
        end
    end
    methods (Test, ParameterCombination='sequential')
        % functions in this method block will only fun if the inputs are
        % all the same length. and each input will be paired
        
    end
end