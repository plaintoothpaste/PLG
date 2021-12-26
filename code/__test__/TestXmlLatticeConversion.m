classdef TestXmlLatticeConversion < matlab.unittest.TestCase
    % tests to check that conversion from xml to lattice and back again is
    % working correctly.
    properties
        xml_default = '__test__/resources/test_cube.xml';
        lattice_default = '__test__/resources/test_cube.lattice';
        xml_tmp = '__test__/resources/test_cube_tmp.xml';
        lattice_tmp = '__test__/resources/test_cube_tmp.lattice';
    end
    methods (Test,TestTags = {'conversion','unit'})
        function test_xml2lattice(testCase)
            %    test_xml2lattice
            
            % run
            xml2lattice(testCase.xml_default,testCase.lattice_tmp,0.1);
            
            % load and verify
            exp_data = loadLattice(testCase.lattice_default);
            act_data = loadLattice(testCase.lattice_tmp);
            testCase.verifyEqual(act_data.vertices, exp_data.vertices);
            testCase.verifyEqual(act_data.struts, exp_data.struts);
            
            % clean up
            delete(testCase.lattice_tmp);
        end
        function test_lattice2xml(testCase)
            %    test_lattice2xml
            % lattice to xml test
            
            % run
            lattice2xml(testCase.lattice_default,testCase.xml_tmp);
            
            % load and verify
            exp_data = loadXml(testCase.xml_default);
            act_data = loadXml(testCase.xml_tmp);
            testCase.verifyEqual(act_data.vertices, exp_data.vertices);
            testCase.verifyEqual(act_data.struts, exp_data.struts);
            
            % clean up
            delete(testCase.xml_tmp);
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
            addpath(genpath('unitCell')); % add classes
        end
    end
    methods (TestClassTeardown)
        % only runs at the end of all the tests
        function cleanup(testCase)
            % remove the test folder and return to starting location
            rmpath(genpath('unitCell')); % add classes
        end
    end
end
function d = loadLattice(file)
    % load a lattice file for testing only
    data = csvread(file);
    num_nodes=data(1,1);
    num_links=data(2,1);
    vert_loc = 3:num_nodes+2;
    strut_loc = num_nodes+3:num_nodes+num_links+2;
            
    d.vertices = data(vert_loc,1:3);
    d.struts   = data(strut_loc,1:2);
end
function d = loadXml(file)
    % load in a xml for testing
    xmlDocument = xmlread(file);
    xmlStructure = parseChildNodes(xmlDocument);
    
    % convert the xml structure to usefull data
    for inc = 1:length(xmlStructure.Children)
        if strcmp(xmlStructure.Children(inc).Name,'vertices')
            vertLoc = inc;
        elseif strcmp(xmlStructure.Children(inc).Name,'struts')
            strutLoc = inc;
        end
    end
    % connections and vertices
    verts = xmlStructure.Children(vertLoc).Children;
    d.vertices = [];
    for inc = 2:2:length(verts)
        d.vertices = [d.vertices;...
            str2double(verts(inc).Attributes(1).Value),...
            str2double(verts(inc).Attributes(2).Value),...
            str2double(verts(inc).Attributes(3).Value)];
    end
    struts = xmlStructure.Children(strutLoc).Children;
    d.struts = [];
    for inc = 2:2:length(struts)
        d.struts = [d.struts;...
            str2double(struts(inc).Attributes(1).Value),...
            str2double(struts(inc).Attributes(2).Value)];
    end
end
function children = parseChildNodes(node)
    % Recurse over node children.
    children = [];
    if node.hasChildNodes
        childNodes = node.getChildNodes;
        numChildNodes = childNodes.getLength;
        allocCell = cell(1, numChildNodes);

        children = struct(             ...
            'Name', allocCell, 'Attributes', allocCell,    ...
            'Data', allocCell, 'Children', allocCell);

        for count = 1:numChildNodes
            theChild = childNodes.item(count-1);
            children(count) = unitCell.makeStructFromNode(theChild);
        end
    end
end