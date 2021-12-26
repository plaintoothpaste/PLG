function xml2lattice(fileIn,fileOut,diameter)
if ~exist('diameter','var')
    diameter=1;
end
if ~exist('fileOut','var')
    [~,fileOut,~] = fileparts(fileIn);
    fileOut = [fileOut,'.lattice'];
end
%% load the xml file
xmlDocument = xmlread(fileIn);
xmlStructure = parseChildNodes(xmlDocument);

%% load the vertices
xml2verts = xmlStructure.Children(2).Children;
test = cellfun(@(x) strcmp(x,'vertex'),{xml2verts(:).Name});
xml2verts(~test)=[];
numVerts = length(xml2verts);
vertices = zeros(numVerts,3);
for inc = 1:numVerts
    x = str2double(xml2verts(inc).Attributes(1).Value);
    y = str2double(xml2verts(inc).Attributes(2).Value);
    z = str2double(xml2verts(inc).Attributes(3).Value);
    vertices(inc,:) = [x,y,z];
end

%% load the struts
xml2struts = xmlStructure.Children(4).Children;
test = cellfun(@(x) strcmp(x,'strut'),{xml2struts(:).Name});
xml2struts(~test)=[];
numStruts = length(xml2struts);
struts = zeros(numStruts,2);
for inc = 1:numStruts
    s1 = str2double(xml2struts(inc).Attributes(1).Value);
    s2 = str2double(xml2struts(inc).Attributes(2).Value);
    struts(inc,:) = [s1,s2];
end

%% save out the data to a lattice file
dlmwrite(fileOut,[numVerts;numStruts]);
data = [vertices,ones(numVerts,1)*diameter];
dlmwrite(fileOut,data,'-append');
data = [struts,ones(numStruts,1)*diameter];
dlmwrite(fileOut,data,'-append');
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
function nodeStruct = makeStructFromNode(node)
% Create structure of node info.

nodeStruct = struct(                        ...
    'Name', char(node.getNodeName),       ...
    'Attributes', unitCell.parseAttributes(node),  ...
    'Data', '',                              ...
    'Children', unitCell.parseChildNodes(node));

if any(strcmp(methods(node), 'getData'))
    nodeStruct.Data = char(node.getData);
else
    nodeStruct.Data = '';
end
end
function attributes = parseAttributes(node)
% Create attributes structure.

attributes = [];
if node.hasAttributes
    theAttributes = node.getAttributes;
    numAttributes = theAttributes.getLength;
    allocCell = cell(1, numAttributes);
    attributes = struct('Name', allocCell, 'Value', ...
        allocCell);
    
    for count = 1:numAttributes
        attrib = theAttributes.item(count-1);
        attributes(count).Name = char(attrib.getName);
        attributes(count).Value = char(attrib.getValue);
    end
end
end