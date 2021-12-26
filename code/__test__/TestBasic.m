classdef TestBasic < matlab.unittest.TestCase
    %Test loading, translating rotating and scaling a lattice.
    properties
        lattice_default = '__test__/resources/test_cube.lattice';
        two_corners = [...
                -0.5,-0.5,-0.5;...
                0.5,0.5,0.5];
        two_index = [1,8]; % the index to get the equivalent corners
    end
    methods (Test,TestTags = {'PLG','unit'})
        function test_loadLattice(testCase)
            %    test_loadLattice
            % test that a lattice file loads correctly
            obj = PLG(testCase.lattice_default);
            
            % validate
            exp_vert = [-0.5,-0.5,0.5];
            act_vert = obj.vertices(2,:);
            testCase.verifyEqual(act_vert,exp_vert);
            
            exp_strut = [1,5];
            act_strut = obj.struts(2,:);
            testCase.verifyEqual(act_strut,exp_strut);
        end
        function test_loadXml(testCase)
            %    test_loadXml
            % test that loading an xml throws an error
            %    test_loadStl
            % test that loading an stl throws an error
            tmp_file = '__test__/tmp.xml';
            copyfile(testCase.lattice_default,tmp_file);
            f = @() PLG(tmp_file);
            testCase.verifyError(f,'PLG:load');
            delete(tmp_file);
            
            tmp_file = '__test__/tmp.stl';
            copyfile(testCase.lattice_default,tmp_file);
            f = @() PLG(tmp_file);
            testCase.verifyError(f,'PLG:load');
            delete(tmp_file);
        end
        function test_loadCsv(testCase)
            %    test_loadCsv
            % test that a warning is issues when loading a csv file
            tmp_file = '__test__/tmp.csv';
            copyfile(testCase.lattice_default,tmp_file);
            f = @() PLG(tmp_file);
            testCase.verifyWarning(f,'PLG:load');
            delete(tmp_file);
        end
        function test_translation(testCase)
            %    test_translation
            % test that translation works as expected
            obj = PLG(testCase.lattice_default);
            
            obj = translate(obj,10,0,0);
            exp_res = testCase.two_corners; 
            exp_res(:,1) = exp_res(:,1)+10;
            act_res = obj.vertices(testCase.two_index,:);
            testCase.verifyEqual(act_res,exp_res);
            
            obj = translate(obj,0,10,0); % second translate
            exp_res = testCase.two_corners; 
            exp_res(:,1) = exp_res(:,1)+10;
            exp_res(:,2) = exp_res(:,2)+10;
            act_res = obj.vertices(testCase.two_index,:);
            testCase.verifyEqual(act_res,exp_res);
            
            obj = translate(obj,0,0,2*pi); % third translate
            exp_res = testCase.two_corners; 
            exp_res(:,1) = exp_res(:,1)+10;
            exp_res(:,2) = exp_res(:,2)+10;
            exp_res(:,3) = exp_res(:,3)+2*pi;
            act_res = obj.vertices(testCase.two_index,:);
            testCase.verifyEqual(act_res,exp_res);
            
            obj = translate(obj,-10,-10,-2*pi); % back to origin
            exp_res = testCase.two_corners; 
            act_res = obj.vertices(testCase.two_index,:);
            testCase.verifyEqual(act_res,exp_res);
        end
        function test_scale(testCase)
            %    test_scale
            % test that scaling the object is perfect
            obj = PLG(testCase.lattice_default);
            obj = scale(obj,10,1,1);
            exp_res = testCase.two_corners; 
            exp_res(:,1) = exp_res(:,1)*10;
            act_res = obj.vertices(testCase.two_index,:);
            testCase.verifyEqual(act_res,exp_res);
        end
        function test_rot(testCase)
            %    test_rot
            obj = PLG(testCase.lattice_default);
            obj = rotate(obj,0,45,0);
            
            new_x = sqrt(2)/2;
            exp_res = testCase.two_corners; 
            exp_res(:,1) = [-new_x;new_x];
            exp_res(:,3) = [0;0];
            
            act_res = obj.vertices(testCase.two_index,:);
            testCase.verifyEqual(act_res,exp_res,'AbsTol',1e-6);
        end
    end
    methods (TestClassSetup)
        % only runs at the start of all the tests
        function setup(testCase)
            % move back to the main directory and create the test folder
            [~,f,~] = fileparts(pwd);
            if strcmp(f,'__test__')
                cd('..');
            end
        end
    end % TestClassSetup
    methods (TestClassTeardown)
        % only runs at the end of all the tests
        function cleanup(testCase)
            % remove the test folder and return to starting location
        end
    end % TestClassTeardown
end

