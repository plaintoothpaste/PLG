# PLG

A group of functions that enable the programmatic lattice generation (PLG).

These can be built up form basic unit cells as well as transformed scaled and joined.

Author(s): Matthew McMillan, David downing, Martin Leary

TODO

*   [critical] radial code fold in.
*   [critical] demo for radial.
*   [medium] better comments in methods.
*   [critical] [paper] fill out tech stack table.

DONE

*   [critical]test demos in the new format.
*   [critical] custom -> .lattice file
*   [critical] update tests to work in the new format.
*   [critical] other BCC

LAST UPDATE: 08.12.2021

[TOC]

## General use

The general overview of use for the PLG is as follows:

1. Define dimensions and replications at minimum using the **set** method.
2. [optional] define origin and base flat, see properties section for the full list.
3. Define the unit cell using the **defineUnit** method.
4. Apply replications using the **cellReplication** method.
5. [recommended] clean duplicate nodes and struts using **cleanLattice**.
6. [optional] perform geometric transformations.
7. [optional] combine with other PLG objects.
8. save the data using a **save** method.

```matlab
% dia, usX,usY,usZ, repX,repY,repZ must be defined to use this example
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
Be aware of the following:

-   Defining a unit cell always overwrites replication information.
-   cellReplication **must** be performed even if the replication is [1,1,1].
-   cleanLattice merges spatially coincident struts and nodes which may lead to different diameters then expected for custom built organic lattices.

## Methods

The following is a list of methods and a quick overview of their function, it should obvious from the name what it is they do.
- PLG: Initialise a PLG object: 0 inputs is an empty object, 1 input loads an existing custom file.
- set: properties are protected (read-only) and can be adjusted using this method. Properties are set using name value pairs.
- defineUnit: used to define the individual unit cell
- cellReplication: replicates a unit cell based on inputs
- cleanLattice: removes duplicate vertices or faces and orders the vertices in Z
- scale: resize the existing structure
- translate: translates the lattice structure.
- rotate: rotates the lattice structure
- plus:   The plus operator `+` is overloaded. This allows the combination of two PLG objects and is usefull for generating composite lattices.
- plot: displays a rendering of the beam model that represents the lattice.
- save: this method saves data to a custom file.

## Available unit cells
The following is a list of presently available unit cells. Custom unit cells can be created from any custom file see **generating a unit cell** for structure of xml file:

| Name               |                            Image                             | Description                                                  |
| :----------------- | :----------------------------------------------------------: | ------------------------------------------------------------ |
| bcc                |            ![BCC](__readMeResources__\bcc_cr.png)            | Cross from the <br>centre of each face.                      |
| centreCross        |    ![centreCross](__readMeResources__\centreCross_cr.png)    | cross from the centre <br>of each face to the centre         |
| horizontalFaceRods | ![horizontalFaceRods](__readMeResources__\horizontalFaceRods_cr.png) | shape on the two planes, top and bottom,<br>that are perpendicular to the z axis. |
| verticalFaceRods   | ![verticalFaceRods](__readMeResources__\VerticalFaceRods_cr.png) | x shape on all the faces that are <br>NOT perpendicular to the z axis. |
| xRods              |          ![xRods](__readMeResources__\xRods_cr.png)          | struts on all edges of the cube<br>that are parallel to the x axis. |
| yRods              |          ![yRods](__readMeResources__\yRods_cr.png)          | struts on all edges of the cube <br>that are parallel to the y axis. |
| zRods              |          ![zRods](__readMeResources__\zRods_cr.png)          | struts on all edges of the cube <br>that are parallel to the z axis. |
| generalCorner      |       ![generalCorner](__readMeResources__\bcc_cr.png)       | a shape that allows for joining of a unit cell at 90 and 45 degrees about all axis where its scaler value is the size of the unit cell you wish to attach too, See complex example. |

## Properties
properties of the PLG are defined using the set method to ensure that only the correct type can be used.
``` matlab
obj = set(obj,'name',value);
```
* resolution - resolution of the struts - scaler integer
* strutDiameter - strut diameter - scaler float
* sphereAddition - determines whether to add spheres to the structure - logical (true or false)
* sphereResolution - required if sphereAddition is true
* sphereDiameter - required if sphereAddition is true
* baseFlat - flattens the spheres at minimum z height vertices to create a flat base for supports - logical (true or false) - does nothing if sphereAddition is false
* unitSize - specifies the size of a single unit in a lattice - 3x1 vector of floats
* replications - specifies the number of copies of the unit cell - 3x1 vector of integers
* origin - starting location for the centre of the initial unit cell - 3x1 vector of floats - default value is [0,0,0]

## Testing

This section needs to be filled in.

## Running in docker

Though it is not recommended the project can be built from a cloned image and then run in docker. This is a 3 step process

1.   Clone and move to the code_ocean branch of the repository.
2.   setup your matlab license in the root of the reposity. **There must be a `license.lic` in root.**
3.   Use `docker image` to create an image from the file `environment/Dockerfile`

```powershell
# 1. create the docker image from the repository, this may take a while.
docker image build ${path to repo}/environment --tag code_ocean 

# 2. IMPORTANT: copy the license.lic before the next step

# 3. To get the docker to run be carefull line returns may mess with powershell
docker run --rm \
-w /code \
-v ${PWD}/code:/code \
-v ${PWD}/data:/data \
-v ${PWD}/results:/results \
-v ${PWD}/license.lic:/MATLAB/licenses/network.lic \ 
code_ocean
```

# custom2stl

A simple function that takes only a input and output. This is used to produce a stl from a custom format.

```matlab
custom2stl(file_in.custom,file_out.stl);
```

# custom2threemf

A simple function that takes only a input and output. This is used to produce a [3mf](https://3mf.io/specification/) file from a custom format.

```matlab
custom2stl(file_in.custom,file_out.3mf);
```

# plotPLG

A class that enables plotting of custom files. The files can then be coloured based on various conditions

```matlab
plotPLG(file_in.custom);
```



# Demos

The following sections contains demonstration files that can be run in MATLAB. These examples are also run automatically in code ocean.

## Regular BCC lattice

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

## Dual density BCZ lattice

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

## Unit cell made from lattices

This is by far the most complex example and creates a FCC structures from other lattices.

```matlab
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
save3mf(obj,'complexLattice.3mf');
saveLattice(obj,'complexLattice.lattice');
%% save out a version with support
obj = addSupport('complexLattice.lattice',diameter/4,0,10,0.1); % add support to all points above minZ to minz (that require it)
obj = padSupport(obj,0.9,diameter/4,0); % extend/add support by a given length with a given strut and ball dia
obj = set(obj,'sphereResolution',12);
obj = set(obj,'resolution',20);
save3mf(obj,'complexLatticeSupported.3mf');
obj = set(obj,'baseFlat',true);
save3mf(obj,'complexLatticeSupportedFlatBall.3mf');
```

![complexExampleWithSupport](__readMeResources__\complexExample.png)

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

## splitStrut
Enables splitting of a bad custom file where beam do interesect but this is not present in the file. splitStruts will identify these and split the beams in two.

## Generating a custom unit cell
unit cells are stored in a xml file to improve flexibility into the future, and minimise confusion with standard output of the program. There are two functions in `code/unitCell` to convert from xml to custom and back again. Not that the xml file stores not diameter data.

```matlab
% fileIn - path to custom file
% fileOut - Optional, path to xml file. Otherwise uses fileIn to make path.
% name - Optional, name placed in the xml file. Otherwise uses fileIn
custom2xml(fileIn,fileOut,name);

% fileIn - path to xml file
% fileOut - Optional, path to custom file. Otherwise uses fileIn to make path.
% diameter - Optional, diameter to apply to the custom file. Otherwise 1 is applied.
xml2custom(fileIn,fileOut,diameter);
```


All the default unit cells are all centred on [0,0,0] and have a bounding box of 1 unit to facilitate scaling. **General corner does not**. While this is not required it is recommended to maximise general application of your unit cell. A comparison between the xml and custom is presented below.

```xml
<?xml version="1.0" encoding="utf-8"?>
<mesh>
   <vertices>
      <vertex x="-0.5" y="-0.5" z="-0.5"/>
      <vertex x="-0.5" y="-0.5" z="0.5"/>
      <vertex x="-0.5" y="0.5" z="-0.5"/>
      <vertex x="-0.5" y="0.5" z="0.5"/>
      <vertex x="0" y="0" z="0"/>
      <vertex x="0.5" y="-0.5" z="-0.5"/>
      <vertex x="0.5" y="-0.5" z="0.5"/>
      <vertex x="0.5" y="0.5" z="-0.5"/>
      <vertex x="0.5" y="0.5" z="0.5"/>
   </vertices>
   <struts name="bcc" type="beam">
      <strut v1="1" v2="5"/>
      <strut v1="2" v2="5"/>
      <strut v1="3" v2="5"/>
      <strut v1="4" v2="5"/>
      <strut v1="5" v2="6"/>
      <strut v1="5" v2="7"/>
      <strut v1="5" v2="8"/>
      <strut v1="5" v2="9"/>
   </struts>
</mesh>
```
 ```
9
8
-0.5, -0.5, -0.5, 1.25
-0.5, -0.5,  0.5, 1.25
-0.5,  0.5, -0.5, 1.25
-0.5,  0.5,  0.5, 1.25
 0.0,  0.0,  0.0, 1.25
 0.5, -0.5, -0.5, 1.25
 0.5, -0.5,  0.5, 1.25
 0.5,  0.5, -0.5, 1.25
 0.5,  0.5,  0.5, 1.25
1, 5, 1
2, 5, 1
3, 5, 1
4, 5, 1
5, 6, 1
5, 7, 1
5, 8, 1
5, 9, 1
 ```



## plgBatch
This is a tool to make a large number of lattice structures. as its input it requires a directory containing at least one json file. Each json file will correspond to a single output file. The code must then be set according to the procedure below:

1. Starting from top to bottom, set all properties. The properties under **TestParameter** will undergo a full factorial design.
2. Scroll down to method and move the curser to the desired output style:
  * runAllCombinations - runs through every single permutation of everything in **TestParameter**.
  * squareUnitCell - uses only unitSizeX for all unit cell dimensions runs through all other **TestParameter**.
  * squareLattice - uses only unitSizeX and repsX
3. click on run current test in the menu bar (ctrl+enter). Alternatively run the following command:
```results = run(plgBatch,'squareLattice')```
4. Be patient matlab internally calculates the full factorial before begining to generate files this may take a while depending on the number of outputs

### Make your own generation function
Generating your own function may be the most usefull for your use case.
1. create a function with the following format:
```function functionNameHere(obj,desiredParameters)```
  * functionNameHere - can be anything not already used and can not be plgBatch or PLG or any PLG functions you plan to call
  * obj - the first input must be the class object itself. this is used to access any constant properties eg: obj.outputFolder
  * desiredParameters - as seperate inputs place any variables in **TestParameter** that you wish to use. New parameters can also be added.
2. follow above instructions

**note:** depending on the number of inputs there can be a big lag between hitting run and the script generating data.
Therefore it is recommended that you test a single output with a script before placing in this class.
