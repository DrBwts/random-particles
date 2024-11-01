# random-particles
Various random space filling codes

# VolFrac.py
A code initially wrtten to be run within Abaqus 6.14-1. The code can be run without Abaqus if all the applicable sections are commented out.
Fills either a rectilinear or spherical volume with randomly placed spheres. The sphere's size's are currently integers I haven't tried floats.
This isnt a closest packing code under the hood it uses an oct-tree as a data structure which lends itself nicely to a cubic lattice packing. Theorically it could be changed to a close packing by mapping the spatial coordinates to the oct-tree cells in some way? 

Be aware that the theoretical limit for the volume fraction for this type of packing is ~0.5236. Values close to this will take a long time to run & may never stop. Values above this will result in the program looping forever (TODO). In the Random_packed_example.jnl a volume fraction of 0.4 was chosen.

The code did once run in FreeCAD 0.16 not sure if it will run in the newer versions. Again the appropriate snippets of code need either commenting or uncommenting. It also used to work quite chunkily with PyPlot & MatPlotLib but I havent updated the code to run with the recent version of MatPlotlib.

Orgininally written in Python 2.17.12 but shoud now work in 3.7
