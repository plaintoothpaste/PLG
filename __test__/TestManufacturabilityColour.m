classdef TestManufacturabilityColour < matlab.unittest.TestCase
    %TESTMANUFACTURABILITYCOLOUR test that components of this class behave
    %as expected
    
    properties
        inputProc = '__test__/resources/basicProcessMap.csv'
        inputCust = '__test__/resources/basicGeom.custom';
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
            actData = manufacturabilityColour.readTables(inData,2,4,6);
            
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
            testCase.verifyError(@()manufacturabilityColour.readTables(inData,2,4,4),'manufacturabilityColour:inconsistentTables');
            
            inData = {'0.2,0,10,20';...
                '0.0001,g,g,g',;...
                '0.5,c,g,g,';...
                '1,y,y,c,';...
                '2,0,10,20';...
                '0.0001,r,r,y';...
                '0.25,r,r,r';...
                '1,r,r,r'};
            testCase.verifyError(@()manufacturabilityColour.readTables(inData,2,4,4),'manufacturabilityColour:inconsistentTables');
        end
        function testReadProcessMap(testCase)
            %    testReadProcessMap
            % test that the reading of csv process maps is correct
            obj = manufacturabilityColour(testCase.inputCust);
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
            obj = manufacturabilityColour(testCase.inputCust);
            actRes = testLengthIncline(obj,testCase.inputProc);
            
            testCase.verifyEqual(actRes.dia,[0.2;0.5;1;4]);
        end
        function testInterpIndex(testCase)
            %    testInterpIndex
            % test interpretation of index number
            incer = [1;4;7];
            obj = manufacturabilityColour(testCase.inputCust);
            actRes = testInterpIndex(obj,testCase.inputProc,incer);
            
            testCase.verifyEqual(actRes.dia,[0.2;0.5;1;4]);
            testCase.verifyEqual(actRes.alpha,[0,10,20,30,90]);
            testCase.verifyEqual(actRes.colour{3,3,1},'c');
            testCase.verifyEqual(actRes.colour{3,2,3},'c');
            testCase.verifyEqual(actRes.colour{4,4,4},'y');
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
            
            % run expansion of 
        end
    end
    methods (TestClassTeardown)
        function teardown(testCase)
            
        end
    end
end

