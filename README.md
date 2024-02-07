# Elle_ice
Elle package with ice

Steps:
1. Generate elle file: input file that has all required ice crystal data for the simulation
2. Generate control script to run a simulation with multiple processes (adapt from past scripts).

## Generate Elle File
1. Download the svg file from https://doi.pangaea.de/10.1594/PANGAEA.933045?format=html#download
2. Use Illustrator to colorcode the grains
3. Convert .svg to .ppm., use customized python codes in 'svg to elle section' at: https://colab.research.google.com/drive/14TMYTLyAIsBpZ_bvM1uMrQeRWZd6bp6x. Or use GIMP (https://www.gimp.org) to convert svg to ppm.
4. Upload ppm to Sherlock, under the folder of stoll_svgppm
5. open ell.simg. convert the ppm to elle file using ```ppm2elle -i input.ppm```. This will open the showelle window. Save the elle file using save as

### Generate Euler angles
1. Open the generated elle file, find the number of grains.
2. Use the sample_watson section in https://colab.research.google.com/drive/14TMYTLyAIsBpZ_bvM1uMrQeRWZd6bp6x to generate watson distribution for the grains
3. Copy the generated euler angle to elle file under the section of ```EULER_3```
4. Note: the euler angles for unodes (mesh points) are not necessary for the initial elle file. But the watson parameter ```k``` can be approximated through the pole plot from: https://doi.pangaea.de/10.1594/PANGAEA.933049?format=html#download

### Generate unodes
1. unodes are the x-y coordinates for the mesh points. It is stored in elle in a point by point form. Each row is one point, the first column is the ID number of the point, the second and the third columns are the x and y coordinates. They are expanded along y direction. i.e. for y for x.
2. Copy the generated unodes into elle file under the section of ```UNODES```.

### Elle file
* EULER_3: the euler angle in degrees for each grain, defined in FLYNNS.
* LOCATION: x - y coordinates of the bnodes
* UNODES: x - y coordinates of the unodes, i.e. the mesh points
* Temperature: in degree C
* BoundaryWidth: 1e-9, unit ??
* MaxNodeSeparation, MinNodeSeparation: unit m
* Unit length: 3e-3, unit ??
* Pressure: 1, unit ??
* SimpleShearOffset & CumulativeSimpleShear: ??
* Timestep: 1.25e7, unit sec ??

Euler Angles IC:
First, use Watson distribution to generate a set of vectors on a sphere in a cartesian coordinate. Second, convert these vectors into Euler angles following Bunge convention. The ranges of the Euler angles (in order) is [0,360], [0,180], [0,360].


### ppc.in
* eqincr, ictrl: ??
* xlfac0, xlfac1: ??
* additional parameters to estimate dislocation: ??

### FNO in julia
It is easier to run on julia notebook
```
using IJulia
notebook()
```
