classdef TestManufacturable < matlab.unittest.TestCase
    % TestManufacturable test that components of this class behave
    % as expected
    
    properties
        inputProc = '__test__/resources/basicProcessMap.csv'
        inputLattice = '__test__/resources/basicGeom.lattice';
        
        inputExcel = '__test__/resources/resourcesData.xlsx';
        
        sheetProc = 'procMap';
        sheetLattice = 'basicGeom';
    end
    
    methods (Test,TestTags = {'unit'})
        function testReadTable(testCase)
            %    testReadTable
            % test that a single table can be converted to usefull data
            inData = {'0.2,0,10,20,30,90';...
                '0.0001,g,g,g,g,g',;...
                '0.5,c,g,g,g,g';...
                '1,y,y,c,g,g';...
                '2,0,10,20,30,90';...
                '0.0001,r,r,y,g,g';...
                '0.5,r,r,r,c,g';...
                '1,r,r,r,c,g'};
            actData = manufacturablePLG.readTables(inData,2,4,6);
            
            testCase.verifyEqual(actData.dia,[0.2;2]);
            testCase.verifyEqual(actData.span,[1e-4;0.5;1]);
        end
        function testReadTableError(testCase)
            %    testReadTableError
            % check that error checking works well
            inData = {'0.2,0,10,20';...
                '0.0001,g,g,g',;...
                '0.5,c,g,g,';...
                '1,y,y,c,';...
                '2,0,10,30';...
                '0.0001,r,r,y';...
                '0.5,r,r,r';...
                '1,r,r,r'};
            testCase.verifyError(@()manufacturablePLG.readTables(inData,2,4,4),'manufacturablePLG:inconsistentTables');
            
            inData = {'0.2,0,10,20';...
                '0.0001,g,g,g',;...
                '0.5,c,g,g,';...
                '1,y,y,c,';...
                '2,0,10,20';...
                '0.0001,r,r,y';...
                '0.25,r,r,r';...
                '1,r,r,r'};
            testCase.verifyError(@()manufacturablePLG.readTables(inData,2,4,4),'manufacturablePLG:inconsistentTables');
        end
        function testReadProcessMap(testCase)
            %    testReadProcessMap
            % test that the reading of csv process maps is correct
            obj = manufacturablePLG(testCase.inputLattice);
            actRes = readProcMap(obj,testCase.inputProc);
            
            testCase.verifyEqual(actRes.dia,[0.2;0.5;1;4]);
            testCase.verifyEqual(actRes.alpha,[0,10,20,30,90]);
            testCase.verifyEqual(actRes.colour{3,3,1},'c');
            testCase.verifyEqual(actRes.colour{3,2,3},'c');
            testCase.verifyEqual(actRes.colour{4,4,4},'y');
        end
        function testStrutLengthAndIncline(testCase)
            %    testStrutLength
            % test that strutLength and incline are correct
            obj = manufacturablePLG(testCase.inputLattice);
            actRes = testLengthIncline(obj,testCase.inputProc);
            
            testCase.verifyEqual(actRes.dia,[0.1;0.175;0.25;0.325;0.4;0.475;0.55;0.625;0.7;0.775]);
            testCase.verifyEqual(actRes.len,[1;0.9;0.8;0.7;0.6;0.7;0.8;0.9;1;1.1],'AbsTol',1e-5);
        end
        function testInterpIndex(testCase)
            %    testInterpIndex
            % test interpretation of index number
            incer = [1;4;7];
            obj = manufacturablePLG(testCase.inputLattice);
            actRes = testInterpIndex(obj,testCase.inputProc,incer);
            
            testCase.verifyEqual(actRes.inclineI,[1;2;4]);
            testCase.verifyEqual(actRes.diaI,[1;1;2]);
            testCase.verifyEqual(actRes.spanI,[3;2;3]);
        end
    end
    methods (Test,TestTags = {'integration'})
        function testBasicRun(testCase)
            % test that the colours array is succefully generated
            obj = manufacturablePLG(testCase.inputLattice);
            obj = runManufacturability(obj,testCase.inputProc);
            
            actColours = obj.colour([1,4,7],:);
            expColours = [1,1,0;0,1,0;0,1,0];
            testCase.verifyEqual(actColours,expColours);
        end
        function testBasicPlot(testCase)
            % test plotting
            obj = manufacturablePLG(testCase.inputLattice);
            obj = runManufacturability(obj,testCase.inputProc);
            
            % warning visual only validation
            plot(obj);
        end
    end
    methods (TestClassSetup)
        function setup(testCase)
            % run location script and setup dirs
            path = cd();
            [~,name,~] = fileparts(path);
            if strcmp(name,'__test__')
                % only move to the root dir if required
                cd('..');
            end
            
            % run expansion of excel file into csv
            data = xlsread(testCase.inputExcel,testCase.sheetLattice);
            data(:,5:end) = [];
            csvwrite(testCase.inputLattice,data);
            
            [~,~,data] = xlsread(testCase.inputExcel,testCase.sheetProc);
            fid = fopen(testCase.inputProc,'w');
            data = data';
            data = data(:);
            TestManufacturable.write2csv(fid,data,6);
        end
    end
    methods (TestClassTeardown)
        function teardown(testCase)
            delete(testCase.inputLattice);
            delete(testCase.inputProc);
        end
    end
    methods (Static)
        function write2csv(fid,data,newLineFreq)
            % recusively write a string or number to a csv where a new line is placed at freq
            if isempty(data)
                fclose(fid);
                return;
            end
            cData = data{1};
            index = mod(length(data),newLineFreq);
            data(1) = [];
            % commer or new line
            
            if index~=1 && ischar(cData)
                fprintf(fid,'%s, ',cData);
            elseif index~=1 && isnan(cData)
                fprintf(fid,',');
            elseif index~=1 && ~ischar(cData) 
                fprintf(fid,'%d, ',cData);
            elseif index==1 && ischar(cData)
                fprintf(fid,'%s,\n',cData);
            elseif index==1 && isnan(cData)
                fprintf(fid,',\n');
            elseif index==1 && ~ischar(cData)
                fprintf(fid,'%d,\n',cData);
            end
            
            TestManufacturable.write2csv(fid,data,newLineFreq);
        end
    end
end

