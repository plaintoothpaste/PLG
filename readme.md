Author(s): Matthew McMillan, 

LAST UPDATE: 06/09/2019

- [description](#description)
  - [available unit cells](#available-unit-cells)
  - [General use](#general-use)
  - [Methods](#methods)
  - [Properties](#properties)
  - [Examples](#examples)
- [Special use cases](#special-use-cases)
  - [Adding Support](#adding-support)
  - [SubClass splitStrut](#subclass-splitstrut)
  - [xlm->custom and custom->xml](#xlm-custom-and-custom-xml)
  - [plgBatch](#plgbatch)

# description
This program takes a variety of inputs and writes a lattice file in various output formats for use in other programs. The use of the plgBatch file simplifies the creation of a large number of lattice(s). Example implementations is present in the complex example file.

## available unit cells
The following is a list of presently available unit cells, **see generating a unit cell** for structure of xml file:

* bcc - cross from each corner into the centre
* ![BCC](__readMeResources__\bcc_cr.png)
* centreCross - cross from the centre of each face to the centre
* ![centreCross](__readMeResources__\centreCross_cr.png)
* horizontalFaceRods - x shape on the two planes, top and bottom, that are perpendicular to the z axis
* ![horizontalFaceRods](__readMeResources__\horizontalFaceRods_cr.png)
* verticalFaceRods - x shape on all the faces that are NOT perpendicular to the z axis
* ![verticalFaceRods](__readMeResources__\VerticalFaceRods_cr.png)
* xRods - strutures on all edges of the cube that are parallel to the x axis
* ![xRods](__readMeResources__\xRods_cr.png)
* yRods - strutures on all edges of the cube that are parallel to the y axis
* ![yRods](__readMeResources__\yRods_cr.png)
* zRods - strutures on all edges of the cube that are parallel to the z axis
* ![zRods](__readMeResources__\zRods_cr.png)
* generalCorner - a shape that allows for joining of a unit cell at 90 and 45 degrees about all axis where its scaler value is the size of the unit cell you wish to attach too, See complex example
* ![generalCorner](__readMeResources__\bcc_cr.png)

current as of commit: ab0b2ead71e301179d2f236f8394eff1b345ae0c

## General use
The general overview of use for the PLG is as follows:

1. Define dimensions and replications at minimum using the **set** method
2. [optional] define origin and base flat, see properties section for the full list
3. Define the unit cell using the **defineUnit** method
4. Apply replications using the **cellReplication** method
5. [recommended] clean duplicate nodes and struts using **cleanLattice**
6. [optional] perform geometric transformations
7. [optional] combine with other PLG objects 
8. [optional] plot the results with the **plot** method
9. save the data using the **save** method

```matlab
obj = PLG() % no input creates an empty object
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
save(obj); % uses a gui save tool
```
Be aware that defining a unit cell will overwrite replication data.
cellReplication must be performed even if the replication is [1,1,1]
cleanLattice merges coincident struts and nodes which may lead to different diameters then expected for custom built organic lattices.

## Methods
The following is a list of methods and a quick overview of their function, it should obvious from the name what it is they do.
- PLG: Creates a PLG object, 1 input loads a custom file.
- set: used to set properties in the PLG using name value pairs
- defineUnit: used to define the individual unit cell
- cellReplication: replicates a unit cell based on inputs
- cleanLattice: removes duplicate vertices or faces and orders the vertices in Z
- scale: resize the existing structure
- translate: translates the lattice structure.
- rotate: rotates the lattice structure
- plus:   combines an existing PLG object into another plg object usefull for generating complex lattice shapes
- plot: displays a rendering of the beam model that represents the lattice.
- save: this method is in its own method group and each sub save method will be named according to its save out type. Eg saveStl saves out stl format.

## Properties
properties of the PLG are defiend using the set method to ensure that only trhe correct type can be used.
``` matlab
obj = set(obj,'name',value);
```
* resolution - resolution of the struts - scaler integer
* strutDiameter - strut diameter - scaler float
* sphereAddition - determines wheter to add spheres to the structure - logical (true or false)
* sphereResolution - required if sphereAddition is true
* sphereDiameter - required if sphereAddition is true
* baseFlat - flattens the spheres at minimum z height vertices to create a flat base for supports - logical (true or false) - does nothing if sphereAddition is false
* unitSize - specifies the size of a single unit in a lattice - 3x1 vector of floats
* replications - specifies the number of copies of the unit cell - 3x1 vector of integers
* origin starting location for the centre of the initial unit cell - 3x1 vector of floats - default value is [0,0,0]

## Examples
A 3x4x5 BCC lattice with x struts, a 0.3mm strut diameter, 4mm unit cell and 0.5mm ball diameter with its origin moved by 6,7,8 and then saved as a stl (12 facet resolution) and 3mf file with a resolution of 30. See complex example for the use of translation, rotation and plus.

```matlab
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
saveStl(obj,'exampleOut.stl');

obj = set(obj,'resolution',30);
obj = set(obj,'sphereResolution',30);
save3mf(obj,'exampleOut.3mf');
```
![flatBase](__readMeResources__\baseFlat.png)

A 3x3x6 5mm unit cell BCZ lattice will be created with a higher density strut diameter of 0.8mm on the top and bottom two row and no z struts in the centre two rows with a 0.5mm diameter.

```matlab
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

saveStl(obj,'exampleOut.stl');
save3mf(obj,'exampleOut.3mf');
```
![plusExample](__readMeResources__\plusExample.png)

# Special use cases

The following use subclasses or are not commonly implemented.

## Adding Support
This class is a submethod of the PLG and enables the addition of support pins.
To use this code you need a fully defined custom file already. in complex example this class is used
This class is not intended for lattice generation but enables the use of the following methods:

* addSupport - class initiation function generates support to the lowest z height - path to custom file, support strutDiameter, support sphereDiameter, critical incline and search range(as a percentage  z height).
* padSupport - extends supports a defined a distance below the minimum z height - object handle, pad height, support strutDiameter, support sphereDiameter.

```matlab
% save out the example above as a custom file
obj = addSupport('path2customWithoutSupport.custom',diameter/4,0,10,0.1);
obj = padSupport(obj,0.9,diameter/4,0);
obj = set(obj,'sphereResolution',12); % required as a custom file does not specify resolution
obj = set(obj,'resolution',20);
save3mf(obj,'geometryWithSupport.3mf');
```
![complexExampleWithSupport](__readMeResources__\complexExample.png)
![complexExampleWithSupportZoom](__readMeResources__\complexExampleZoom.png)

## SubClass splitStrut
Enables splitting of a bad custom file where beam do interesect but this is not present in the file. splitStruts will identify these and split the beams in two.

## xlm->custom and custom->xml
unit cells are stored in a xml file to improve flexibility into the future and minimise confusion with standard output of the program. There are functions to convert from xml to custom and back again.
Current unit cells are all centred on [0,0,0] and have a bounding box of 1 unit to facilitate scaling. **General corner does not**. While this is not required it is recommended to maximise general application.

 ![xml2custom](__readMeResources__\customVsXml.PNG)

## plgBatch
1. starting from top to bottom set all properties. The properties under **TestParameter** will undergo a full factorial design.
2. Scroll down to method and move the curser to the desired output style:
  * runAllCombinations - runs through every single permutation of everything in **TestParameter**.
  * squareUnitCell - uses only unitSizeX for all unit cell dimensions runs through all other **TestParameter**.
  * squareLattice - uses only unitSizeX and repsX
3. click on run current test in the menu bar (ctrl+enter). Alternatively run the following command:
```results = run(plgBatch,'squareLattice')```
4. Be patient matlab internally calculates the full factorial before begining to generate files this may take a while depending on the number of outputs

### make your own generation function
Generating your own function may be the most usefull for your use case.
1. create a function with the following format:
```function functionNameHere(obj,desiredParameters)```
  * functionNameHere - can be anything not already used and can not be plgBatch or PLG or any PLG functions you plan to call
  * obj - the first input must be the class object itself. this is used to access any constant properties eg: obj.outputFolder
  * desiredParameters - as seperate inputs place any variables in **TestParameter** that you wish to use. New parameters can also be added.
2. follow above instructions

**note:** depending on the number of inputs there can be a big lag between hitting run and the script generating data.
Therefore it is recommended that you test a single output with a script before placing in this class.
