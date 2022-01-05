%Create an image of all unit cells in the unit cell folder as well as any
% additional combination unit cells listed below. These images are also 
% used in the readme.

% names ralative to the unit cell path
unit_names = {
    {'verticalFaceRods','horizontalFaceRods','xRods','yRods','zRods'};
    {'xRods','yRods','zRods'};
    };

[code_location, unit_cell_location, output_location] = pathSetup();
addpath(code_location);

unit_names = getAllUnitCells(unit_cell_location,unit_names);
for i = 1:length(unit_names)
    unit = unit_names{i};
    genericPlot(unit,output_location)
end
rmpath(code_location);


function [code_location, unit_cell_location, output_location] = pathSetup()
    %    pathSetup
    % this code can be run from mulitple locations so setup paths for 
    % adding directories etc.
    pathPLG = which('plotUnitCell');
    runLocation = fileparts(pathPLG);

    code_location = [runLocation , '/../code'];
    unit_cell_location = [runLocation , '/../code/unitCell'];
    output_location = runLocation;
end

function unit_names = getAllUnitCells(unit_cell_location,unit_names)
    % Add to all the existing xml unit cells
    dir_xml = dir([unit_cell_location , '/*.xml']);
    
    for i = 1:length(dir_xml)
        n = split(dir_xml(i).name,'.');
        unit_names{end+1} = {n{1}};
    end
end

function genericPlot(unit_cell,out_directory)
    % A generic thumbnail plot of a given unit cell.
    %    The file name is auto generated based on the inputs.
    obj=PLG();
    obj=set(obj,'strutDiameter',1);
    obj=set(obj,'sphereDiameter',1);
    obj=set(obj,'unitSize',[1,1,1]);
    obj=set(obj,'replications',[1,1,1]);
    obj=defineUnit(obj,unit_cell);
    obj=cellReplication(obj);
    
    [f,a] = plot(obj);
    view([-5 2 1]);
    f.Units = 'centimeters';
    f.Position=[3.5, 3.5, 4, 4];
    f.PaperPositionMode = 'auto';
    a.FontSize = 7;
    grid("on");
    h=get(a,'Children');
    set(h,'LineWidth',2);

    % generate name and save
    file_name = [out_directory,filesep];
    for i = 1:length(unit_cell)
        file_name = [file_name,unit_cell{i},'_'];
    end
    file_name = [file_name(1:(end-1)),'.png'];
    print(f,file_name,'-dpng','-r300');

    close(f);
end

