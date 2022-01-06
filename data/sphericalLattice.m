% build a custom spherical lattice structure for the SLM process.
% see https://bitbucket.org/plg-lattice/spherical_lattice for full version
function sphericalLattice()
    % create a sperical lattice
    
    %% paths and configuration
    %there are a large number of settings so they are stored in a json file
    lattice_path = fileparts(which('sphericalLattice'));
    spec = jsondecode(fileread([lattice_path,filesep,'spec.json']));
    % add paths to json structure
    spec.general.path.lattice = lattice_path; 
    spec.general.path.PLG = [spec.general.path.lattice,'/../code'];
    spec.general.path.saveDirectory = [spec.general.path.lattice,'/../results'];
    spec.general.path.processMap = [spec.general.path.lattice,'/processMap.csv'];
    addpath(spec.general.path.PLG);
    
    %% build segments
    % each part of the sphere is made up of segments build them in isolation
    [segment_names,segment_data] = segmentRetrieve(spec);
    
    % only required if rebuilding segment, see json for setting
    for inc = 1:length(segment_names)
        buildSection(segment_names{inc},segment_data{inc},spec.general);
        makeImage(segment_names{inc},spec.general);
    end
    
    %% support
    % Add support pins.
    buildSupport(spec.support,spec.general);
    makeImage(spec.support.name,spec.general);
    
    %% full build section
    % put the above sections together
    name_list = [segment_names;spec.support.name];
    combineSections(name_list,spec.file_name.segment,spec.file_name.full,spec.general);

    % plot manufacturability
    makeImage(spec.file_name.segment,spec.general);
    makeImage(spec.file_name.full,spec.general);

    %calculate statistics
    makeStatistics(spec.file_name.full,spec.general);

    % make stl file
    makeStl(spec.file_name.segment,spec.file_name.stl_segment,spec.general);
    makeStl(spec.file_name.full,spec.file_name.stl_full,spec.general);

    rmpath(spec.general.path.PLG);
end


%% functions
function [name,segment] = segmentRetrieve(spec)
    % from the global json data retrieve the list of file names and segments
    
    c_name = fieldnames(spec.segments);
    if isempty(c_name)
        name = [];
        segment = [];
        return;
    end
    c_name = c_name{1};
    c_segment = spec.segments.(c_name);
    spec.segments = rmfield(spec.segments,c_name);
    c_name = {sprintf('%s_%d_%d.lattice',c_name,c_segment.seg(1),c_segment.seg(2))};
    
    % recurse
    [f_name,f_segment] = segmentRetrieve(spec);
    name = [c_name;f_name];
    segment = [{c_segment},f_segment];
end % segmentRetrieve

function [name,core] = coreRetrieve(names,spec)
    % from the global json data retrieve the list of file names and segments
    
    if isempty(names)
        name = [];
        core = [];
        return;
    end
    c_name = names{1};
    names(1) = [];
    segment_names = fieldnames(spec.segments);
    
    c_core = spec.core.(c_name);
    test = strcmp(c_core.derived,segment_names);
    segment = spec.segments.(segment_names{test});
    c_core.seg = segment.seg;
    c_core.reps = segment.reps;
    c_core.mirror = segment.mirror;
    spec.core = rmfield(spec.core,c_name);
    c_name = {sprintf('%s_%d_%d.lattice',c_name,c_core.seg(1),c_core.seg(2))};
    
    % recurse
    [f_name,f_core] = coreRetrieve(names,spec);
    name = [c_name;f_name];
    core = [{c_core},f_core];
end % coreRetrieve

function makeImage(f_name,g)
    % saves a figure at a standard view and size
    obj = manufacturablePLG([g.path.saveDirectory,filesep,f_name]);
    obj = runManufacturability(obj,g.path.processMap);
    [f,a] = plot(obj);
    
    [~,f_name] = fileparts(f_name);
    a.View = [0,0];
    f.Units = 'centimeters';
    f.Position=[3.5, 3.5, 12, 12];
    f.PaperPositionMode = 'auto';
    a.FontSize = 7;
    h=get(a,'Children');
    set(h,'LineWidth',2);
    set(h,'MarkerSize',2);
    grid("on");

    print(f,[g.path.saveDirectory,filesep,f_name,'.png'],'-dpng','-r300');
    close(f);
end % makeImage

function obj = buildSection(f_name,spec,g)
    % generate a portion of the sperical lattice
    obj = sphericalPLG();
    obj = set(obj,'strutDiameter',g.diameter);
    obj = set(obj,'sphereDiameter',g.diameter);
    obj = set(obj,'sphereAddition',true);
    r = spec.seg*pi/180; % convert to radians
    obj = set(obj,'unitSize',[g.unit_size,2*pi/(g.reps.tangential),(r(2)-r(1))/spec.reps]);
    obj = set(obj,'replications',[g.reps.radial,1,spec.reps]);
    
    % unit cell
    [~,f,ext] = fileparts(spec.unit{1});
    if ~isempty(ext) || strcmp(ext,'.lattice')
        p = [g.path.PLG,'/unitCell'];
        addpath(p);
        pp = [p,'/',f,'.xml'];
        lattice2xml([g.path.lattice,filesep,spec.unit{1}],pp);
        rmpath(p);
        [~,s,~] = fileparts(pp);
        spec.unit = {s};
    end
    obj = defineUnit(obj,spec.unit);
    % replicate and translate
    obj = cellReplication(obj); % requires the replications to be set
    obj = translate(obj,g.R + g.unit_size*(1/2-g.reps.radial), pi/g.reps.tangential,r(1) + 1/2*(r(2)-r(1))/spec.reps-pi/2);
    obj = cart2polar(obj);
    
    % duplicate the cap
    if spec.mirror
       obj2 = rotate(obj,180,0,0);
       obj = obj+obj2;
    end
    obj = cleanLattice(obj,g.tolerance); % remove duplicate points and merge close points
    saveLattice(obj,[g.path.saveDirectory,filesep,f_name]);
    
    if ~isempty(ext) || strcmp(ext,'.lattice')
        delete(pp)
    end
end % buildSection

function buildSupport(sup,g)
    % build some support pins
    
    % combine the files which need to be supported
    f_namer = @(p,s) sprintf('%s/%s*.lattice',p,s);
    n = length(sup.parts);

    s = dir(f_namer(g.path.saveDirectory,sup.parts{1}));
    obj = PLG([s.folder,filesep,s.name]);
    if n>1
        for inc = 2:n
            s = dir(f_namer(g.path.saveDirectory,sup.parts{inc}));
            obj = obj + PLG([s.folder,filesep,s.name]);
        end
    end
    obj = cleanLattice(obj,g.tolerance);
    combined_support_name = [g.path.saveDirectory,filesep,sup.name];
    saveLattice(obj,combined_support_name);
    
    % load the combined file and add supports
    obj = addSupport(combined_support_name,sup.diameter,sup.diameter,sup.incline,sup.search);
    obj = padSupport(obj,sup.pad,sup.diameter,sup.diameter);
    obj = cleanLattice(obj);
    
    saveLattice(obj,combined_support_name);
end % buildSupport

function combineSections(file_names,file_out_seg,file_out_full,g)
    % combine multiple segments and revolve to generate the entire sphere.
    n = length(file_names);
    output_dir = [g.path.saveDirectory,filesep];
    obj = PLG([output_dir,file_names{1}]);
    for inc = 2:n
        obj = obj + PLG([output_dir,file_names{inc}]);
    end
    obj = cleanLattice(obj,g.tolerance);
    saveLattice(obj,[output_dir,file_out_seg]);
    
    % revolve to get whole sphere
    rot = 360/g.reps.tangential;
    new = obj;
    all = obj;
    for inc = 1:g.reps.tangential
        new = rotate(new,0,0,rot);
        all = all + new;
    end
    all = cleanLattice(all,g.tolerance);
    saveLattice(all,[output_dir,file_out_full]);
end % combineSections

function makeStl(file_name_in,file_name_out,g)
    output_dir = [g.path.saveDirectory,filesep];
    obj = PLG([output_dir,file_name_in]);
    obj = set(obj,'resolution',g.resolution);
    obj = set(obj,'sphereResolution',g.resolution);
    obj = set(obj,'sphereAddition',true);
    obj = set(obj,'baseFlat',true);
    saveStl(obj,[output_dir,file_name_out]);
end % makeStl

function makeStatistics(f_name,g)
    % save out a statistics excel based on f_name
    obj = statisticsPLG([g.path.saveDirectory,filesep,f_name],g.path.processMap);
    
    [~,f_name] = fileparts(f_name);
    save(obj,[g.path.saveDirectory,filesep,f_name,'.xlsx']);
end