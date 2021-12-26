function lattice2xml(fileIn,fileOut,name)
% takes a lattice input file and converts it the the xml storage file used
% by the PLG
if ~exist('fileOut','var')
    [~,fileOut,~] = fileparts(fileIn);
    fileOut = [fileOut,'.xml'];
end
if ~exist('name','var')
    [~,name,~] = fileparts(fileIn);
end

%% load the lattice file
data = csvread(fileIn);
numNodes=data(1,1);
numLinks=data(2,1);
vertices = data(3:numNodes+2,1:3);
struts   = data(numNodes+3:numNodes+numLinks+2,1:2);
clear data
%% set up the xml file
xmlDocument = com.mathworks.xml.XMLUtils.createDocument('mesh');
mesh_xml = xmlDocument.getDocumentElement;
vertices_xml = xmlDocument.createElement('vertices');
mesh_xml.appendChild(vertices_xml);

struts_xml = xmlDocument.createElement('struts');
struts_xml.setAttribute('name',name);
struts_xml.setAttribute('type','beam');
mesh_xml.appendChild(struts_xml);

%% write the vertices to the xml file
for inc = 1:numNodes
    currentStrut = vertices(inc,:);
    vertex = xmlDocument.createElement('vertex');
    vertex.setAttribute('x',num2str(currentStrut(1)));
    vertex.setAttribute('y',num2str(currentStrut(2)));
    vertex.setAttribute('z',num2str(currentStrut(3)));
    vertices_xml.appendChild(vertex);
end
%% write the struts to the xml file
for inc = 1:numLinks
    currentStrut = struts(inc,:);
    strut = xmlDocument.createElement('strut');
    strut.setAttribute('v1',num2str(currentStrut(1)));
    strut.setAttribute('v2',num2str(currentStrut(2)));
    struts_xml.appendChild(strut);
end

%% save out the xml file
xmlwrite(fileOut,xmlDocument);

end