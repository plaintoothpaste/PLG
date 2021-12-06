function results = runTests(tag,verbose)
%runTests creates an ojbect to run tests.
%    Runs all the tests for PLG within matlab and return the results with options for tagging and verbosity.
% Inputs
%    tag    : Optional for filtering tests.
%    verbose: Optional boolean that if set to true displays more information for tests.
    if ~exist('tag','var')
        tag = '*'; % all
    end
    if ~exist('verbose','var')
        verbose = true;
    end
    import matlab.unittest.TestSuite;
    import matlab.unittest.TestRunner;
    if verbose
        runner = TestRunner.withTextOutput;
    else
        runner = TestRunner.withNoPlugins;
    end
    suite = matlab.unittest.TestSuite.fromFolder('__test__','Tag',tag); % run all tags
    results = run(runner,suite);
    results = table(results);
end