function autoTestRunner(tag)
% automatically run tests based on tag input
% if no input the all tests are run

import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.TAPPlugin;
import matlab.unittest.plugins.ToFile;

% if no tag supplied
if ~exist('tag','var')
    tag = '*';
end

% remove old tap file if present
tapFileLocation = fullfile('__test__/resources/testResults.tap');
if ~isempty(dir(tapFileLocation))
    delete(tapFileLocation);
end

% setup suite and randomise the order to ensure there is no interdependence
suite = matlab.unittest.TestSuite.fromFolder('__test__','Tag',tag);
suite = suite(randperm(numel(suite))); 

runner = TestRunner.withTextOutput();
tapFile = tapFileLocation;
runner.addPlugin(TAPPlugin.producingOriginalFormat(ToFile(tapFile)));

results = runner.run(suite);
display(results);
end