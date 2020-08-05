function results = runTests(tag,verbose)
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