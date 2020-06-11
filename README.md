# random-particles
Various random space filling codes

# VolFrac.py
A code initially wrtten to be run within Abaqua 6.14-1. The code can be run without Abaqus if all the applicable sections are commented out.
Fills either a rectilinear or spherical volume with randomly placed spheres. The sphere's size's are currently integers haven't tried floats.
This sint a closest packing code under the hood it uses an oct-tree as a data structure which lends itself nicely to a cubic lattice packing/ Theorically it could be changed to a close packing by mapping the spatial coordinates to the oct-tree cells in some way?

The code did once run in FreeCAD 0.16 not sure if it run in the newer versions. Again the appropriate snippets of code need either commenting or uncommenting. It also used to work quite chunkily with PyPlot & MatPlotLib but I havent updated the code to run with the recent version of MatPlotlib.

Orgininally written in Python 2.17.12 but shoud now work in 3.7
