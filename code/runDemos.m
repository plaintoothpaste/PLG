%% A script that runs all demo files in the readme
mostBasicDemo();
regularBccLattice();
dualDensityBccLattice();
complexExample();
addSupportToComplexExample();
radialLattice();
sphereLattice(); % warning this example is larger then all the others put together.

%% each demo is in its own function below
function mostBasicDemo()
    dia = 1;
    us = 10;
    usX = us; usY = us; usZ = us;
    reps = 2;
    repX = reps; repY = reps; repZ = reps; 
    
    obj = PLG(); % no input creates an empty object
    obj = set(obj,'resolution',20);
    obj = set(obj,'sphereResolution',12);
    obj = set(obj,'strutDiameter',dia);
    obj = set(obj,'sphereDiameter',1.5*dia);
    obj = set(obj,'sphereAddition',true); % default value is false
    obj = set(obj,'unitSize',[usX,usY,usZ]);
    obj = set(obj,'replications',[repX,repY,repZ]);
    obj = defineUnit(obj,{'bcc','zRods'});
    obj = cellReplication(obj); % requires the replications to be set
    obj = cleanLattice(obj); % remove any vertices and struts that are coincident
    % save(obj); % uses a gui save tool, this is removed in this version for automation reasons
    save3mf(obj,'../results/mostBasicDemo.3mf');

    plotDemo(obj,'../results/mostBasicDemo.png');
end

function regularBccLattice()
    obj = PLG();
    obj = set(obj,'resolution',12);
    obj = set(obj,'strutDiameter',0.3);
    obj = set(obj,'unitSize',[4,4,4]);
    obj = set(obj,'sphereAddition',true);
    obj = set(obj,'baseFlat',true);
    obj = set(obj,'sphereResolution',12);
    obj = set(obj,'sphereDiameter',0.5);
    obj = set(obj,'origin',[6,7,8]);
    obj = set(obj,'replications',[3,4,5]);
    obj = defineUnit(obj,{'bcc','xRods'});
    obj = cellReplication(obj);
    obj = cleanLattice(obj);
    saveStl(obj,'../results/regularBccLattice.stl');
    
    obj = set(obj,'resolution',30);
    obj = set(obj,'sphereResolution',30);
    save3mf(obj,'../results/regularBccLattice.3mf');

    plotDemo(obj,'../results/regularBccLattice.png');
end

function dualDensityBccLattice()
    bottomLayer = PLG();
    bottomLayer = set(bottomLayer,'resolution',12);
    bottomLayer = set(bottomLayer,'strutDiameter',0.8);
    bottomLayer = set(bottomLayer,'unitSize',[5,5,5]);
    bottomLayer = set(bottomLayer,'sphereAddition',true);
    bottomLayer = set(bottomLayer,'sphereResolution',12);
    bottomLayer = set(bottomLayer,'sphereDiameter',0.8);
    bottomLayer = set(bottomLayer,'origin',[0,0,0]);
    bottomLayer = set(bottomLayer,'replications',[3,3,2]);
    
    topLayer = bottomLayer;
    topLayer = set(topLayer,'origin',[0,0,20]);
    
    midLayer = bottomLayer;
    midLayer = set(midLayer,'strutDiameter',0.5);
    midLayer = set(midLayer,'sphereDiameter',0.5);
    midLayer = set(midLayer,'origin',[0,0,10]);
    
    bottomLayer = defineUnit(bottomLayer,{'bcc','zRods'});
    bottomLayer = cellReplication(bottomLayer);
    bottomLayer = cleanLattice(bottomLayer);
    topLayer = defineUnit(topLayer,{'bcc','zRods'});
    topLayer = cellReplication(topLayer);
    topLayer = cleanLattice(topLayer);
    midLayer = defineUnit(midLayer,{'bcc'});
    midLayer = cellReplication(midLayer);
    midLayer = cleanLattice(midLayer);
    
    obj = bottomLayer + topLayer + midLayer;
    obj = cleanLattice(obj);
    
    saveStl(obj,'../results/dualDensityBccLattice.stl');
    save3mf(obj,'../results/dualDensityBccLattice.3mf');

    plotDemo(obj,'../results/dualDensityBccLattice.png');
end

function complexExample()
    unitSize = 4;
    repSpacing = 160;
    cornerSpacing = repSpacing+2*(0.5*unitSize+sqrt(unitSize^2/2));
    diameter = 0.5;
    originOffset = unitSize+sqrt(unitSize^2/2);
    %% make balls
    obj = PLG;
    obj = set(obj,'resolution',20);
    obj = set(obj,'strutDiameter',diameter);
    obj = set(obj,'sphereAddition',true);
    obj = set(obj,'sphereDiameter',1.5*diameter);
    obj = set(obj,'sphereResolution',12);
    obj = set(obj,'unitSize',[unitSize,unitSize,unitSize]);
    obj = set(obj,'origin',[0,0,0]);
    obj = defineUnit(obj,{'generalCorner'});
    obj = set(obj,'replications',[1,1,1]);
    obj = cellReplication(obj);
    obj = cleanLattice(obj);
    
    obj1 = translate(obj,0,0,cornerSpacing);
    obj2 = translate(obj,cornerSpacing/2,0,cornerSpacing/2);
    obj3 = translate(obj,cornerSpacing,0,0);
    obj4 = translate(obj,cornerSpacing,0,cornerSpacing);
    obj = obj1+obj2+obj3+obj4+obj; % 1 face
    
    obj1 = rotate(obj,0,0,90);
    obj2 = translate(obj,0,cornerSpacing,0);
    obj3 = translate(obj1,cornerSpacing,0,0);
    obj = obj1+obj2+obj3+obj;
    cornerObj = cleanLattice(obj);
    %% make the straight sections
    %along x axis
    obj = PLG;
    obj = set(obj,'resolution',20);
    obj = set(obj,'strutDiameter',diameter);
    obj = set(obj,'sphereAddition',true);
    obj = set(obj,'sphereDiameter',diameter);
    obj = set(obj,'sphereResolution',12);
    obj = set(obj,'unitSize',[unitSize,unitSize,unitSize]);
    obj = set(obj,'origin',[originOffset,0,0]);
    obj = defineUnit(obj,{'bcc','xRods'});
    obj = set(obj,'replications',[repSpacing/unitSize,1,1]);
    obj = cellReplication(obj);
    
    obj1 = translate(obj,0,0,cornerSpacing);
    obj = obj1+obj;
    obj1 = translate(obj,0,cornerSpacing,0);
    objX = obj1+obj;
    
    % along z axis, rotate x by 90
    objY = rotate(objX,0,0,90);
    objY = translate(objY,cornerSpacing,0,0);
    
    % along y axis, rotate x by 90
    objZ = rotate(objX,0,90,0);
    objZ = translate(objZ,0,0,cornerSpacing);
    %plot(objX+objY+objZ);

    %% make diagonals
    obj = PLG;
    obj = set(obj,'resolution',20);
    obj = set(obj,'strutDiameter',diameter);
    obj = set(obj,'sphereAddition',true);
    obj = set(obj,'sphereDiameter',diameter);
    obj = set(obj,'sphereResolution',12);
    specUnitSize = (sqrt(2*(repSpacing/2)^2)-unitSize*2)/20;
    obj = set(obj,'unitSize',[specUnitSize,unitSize,unitSize]);
    specialOffset = specUnitSize/2+unitSize/2+sqrt(unitSize^2/2);
    obj = set(obj,'origin',[specialOffset,0,0]);
    obj = defineUnit(obj,{'bcc','xRods'});
    obj = set(obj,'replications',[21,1,1]);
    obj = cellReplication(obj);
    obj = cleanLattice(obj);
    obj1 = rotate(obj,0,-45,0);
    obj2 = rotate(obj,0,-45+90,0);
    obj3 = rotate(obj,0,-45+180,0);
    obj4 = rotate(obj,0,-45+270,0);
    obj = translate(obj1+obj2+obj3+obj4,cornerSpacing/2,0,cornerSpacing/2);
    
    diagXZ = obj + translate(obj,0,cornerSpacing,0);
    diagYZ = rotate(diagXZ,0,0,90);
    diagYZ = translate(diagYZ,cornerSpacing,0,0);
    % plot(diagXZ+diagYZ);
    
    %% combine all the objects and save out
    obj = objX+objY+objZ+cornerObj+diagXZ+diagYZ;
    obj = cleanLattice(obj);
    obj = set(obj,'sphereResolution',20);
    save3mf(obj,'../results/complexLattice.3mf');
    saveLattice(obj,'../results/complexLattice.lattice');

    plotDemo(obj,'../results/complexLattice.png');
end

function addSupportToComplexExample()
    % save out a version of the complex example with support
    %    requires complexExample to be run
    diameter = 0.5;
    obj = addSupport('../results/complexLattice.lattice',diameter/4,0,10,0.1); % add support to all points above minZ to minz (that require it)
    obj = padSupport(obj,0.9,diameter/4,0); % extend/add support by a given length with a given strut and ball dia
    obj = set(obj,'sphereResolution',12);
    obj = set(obj,'resolution',20);
    save3mf(obj,'../results/complexLatticeSupported.3mf');
    obj = set(obj,'baseFlat',true);
    save3mf(obj,'../results/complexLatticeSupportedFlatBall.3mf');

    plotDemo(obj,'../results/complexLatticeSupportedFlatBall.png');
end

function radialLattice()
    % create a lattice in a radial format
    obj = radialPLG(); % no input creates an empty object
    obj = set(obj,'resolution',6);
    obj = set(obj,'sphereResolution',6);
    obj = set(obj,'strutDiameter',1);
    obj = set(obj,'sphereDiameter',1.25);
    obj = set(obj,'sphereAddition',true); % default value is false
    obj = set(obj,'unitSize',[2,2*pi/12,2]);
    obj = set(obj,'replications',[2,12,2]);
    obj = defineUnit(obj,{'bcc','zRods'});
    obj = cellReplication(obj); % requires the replications to be set
    obj = translate(obj,12,0,0);
    obj = cart2radial(obj);
    saveStl(obj,'../results/radial.stl');

    plotDemo(obj,'../results/radial.png');
end

function sphereLattice()
    addpath('../data');

    sphericalLattice();

    rmpath('../data');
end

function plotDemo(obj, file_name)
    % plot a image of a demo file
    [f,a] = plot(obj);
    view([-5, 2, 1]);
    f.Units = 'centimeters';
    f.Position=[3.5, 3.5, 6, 6];
    f.PaperPositionMode = 'auto';
    a.FontSize = 7;
    grid("on");
    h=get(a,'Children');
    set(h,'LineWidth',2);
    set(h,'MarkerSize',2);
    
    print(f,file_name,'-dpng','-r300');
    close(f);
end