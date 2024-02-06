# Elle_ice
Elle package with ice

Steps:
1. Generate elle file: input file that has all required ice crystal data for the simulation
2. Generate shelle file: control script to run a simulation with multiple processes.

## Generate Elle File
1. Download the svg file from https://doi.pangaea.de/10.1594/PANGAEA.933045?format=html#download
2. Use Illustrator to colorcode the grains
3. Convert .svg to .ppm., use customized python codes in svg to elle section at: https://colab.research.google.com/drive/14TMYTLyAIsBpZ_bvM1uMrQeRWZd6bp6x
4. Upload ppm to Sherlock, under the folder of stoll_svgppm
5. open ell.simg. convert the ppm to elle file using '''ppm2elle -i input.ppm'''

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
