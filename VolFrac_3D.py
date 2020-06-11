# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 15:42:02 2015

@author: Nicholas Roberts
"""

import numpy as np
import copy
from collections import OrderedDict

# Uncomment if using matplotlib & pyplot 
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

# Abaqus imports, comment out if not running in Abaqus
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

# FreeCAD imports, comment out if not running FreeCAD
#import Part, FreeCAD, math, copy
#from FreeCAD import Base

class VolumeFill():
    
    def __init__(self, dims, diameters, isSphere=False):
        
        self.depth = len(diameters)
        self.diameters = diameters
        self.dims = dims
        self.majorRadius = np.max(dims)/2.0
        self.isSphere = isSphere
        self.allCentroids = []
        self.allCoordinates = []
        self.tree = self.treeSetUp()
        
        
    def treeSetUp(self):       
        '''Set tree's initial structure
        For each level in treeArray
        [[node number, full, leaf, parent, [children], coordinates, [centroid]],   [node number, full, leaf, parent, [children], coordinates, [centroid]], ...]'''
        
        # Create root level & node    
        leaf = False    
        numChildren = (self.dims[0] * self.dims[1] * self.dims[2])/self.diameters[0]**3
        children = [(kids) for kids in range(numChildren)]
        rootCoords = np.array([[0, 0, 0], [0, self.dims[0], 0], [self.dims[0], 0, 0], [self.dims[0], self.dims[1], 0], [0, 0, self.dims[2]], [0, self.dims[1], self.dims[2]], [self.dims[0], 0, self.dims[2]], [self.dims[0], self.dims[1], self.dims[2]]])
        xCentroid = np.sum(rootCoords[:,0])/8.0
        yCentroid = np.sum(rootCoords[:,1])/8.0
        zCentroid = np.sum(rootCoords[:,2])/8.0
        rootCentroid = [xCentroid, yCentroid, zCentroid]
        rootNode = [0, False, copy.copy(leaf), None, copy.copy(children), rootCoords, rootCentroid]
        
        treeArray = []
        treeArray.append(rootNode)
        self.allCoordinates.append(rootCoords)
        self.allCentroids.append(rootCentroid)
        
        
        # Create level 1 nodes
        level = []
        levelCentroids, levelCoordinates = self.coordGenerator(1, self.dims, [0,0,0])
        self.allCentroids.append(copy.copy(levelCentroids))
        self.allCoordinates.append(levelCoordinates)
        
        for child in treeArray[0][4]:
        
            if self.depth == 1:
                leaf = True

            childNode = [child, False, leaf, 0, [], levelCoordinates[child], levelCentroids[child]]
            level.append(copy.deepcopy(childNode))
            
        treeArray.append(copy.deepcopy(level))
        
        
        # Create other level nodes
        for treeLevel in range(1, self.depth):
            
            level = []
            childCount = 0
            
            for cInx, child in enumerate(treeArray[treeLevel]):
                
                if treeLevel == self.depth - 1:
                    leaf = True
                    
                childDims = [int(np.max(child[5][:,0])-np.min(child[5][:,0])), int(np.max(child[5][:,1])-np.min(child[5][:,1])), int(np.max(child[5][:,2])-np.min(child[5][:,2]))]
                childOffset = [int(np.min(child[5][:,0])), int(np.min(child[5][:,1])), int(np.min(child[5][:,2]))]
                childCentroids, childCoords = self.coordGenerator(treeLevel+1, childDims, childOffset)
                
                for childCell in zip(childCoords, childCentroids):
                    
                    distanceFromCentre = self.euclidean(rootCentroid, childCell[1])
                    if self.isSphere and distanceFromCentre > self.majorRadius and treeLevel == self.depth - 1:
                        full = True
                    else:
                        full = False
                        
                    childNode = [childCount, full, leaf, child[0], [], childCell[0], childCell[1]]
                    level.append(copy.deepcopy(childNode))
                    treeArray[treeLevel][cInx][4].append(childCount)
                    childCount += 1
#                    print treeLevel, childCount
                    
            treeArray.append(copy.deepcopy(level))
            
        
        return treeArray
        
    def coordGenerator(self, level, dims, offset):     
        '''Generate the coordinates for the vertices of every child cell within a cell'''
        
        # Generate 2D mesh intervals
        levelCoords = np.array([(x+offset[0], y+offset[1]) for x in range(0, dims[0]+1, self.diameters[level-1]) for y in range(0, dims[1]+1, self.diameters[level-1])]) 
        gap = int(np.sqrt(len(levelCoords)))
        intervals = np.array([(temp) for temp in range(0, len(levelCoords)-gap, gap)])  
        
        # Generate 2D mesh
        cells2Dsimplices = []
        for idx in intervals:
            
            tempSimp = np.reshape(np.array([levelCoords[idx:idx+2,:], levelCoords[idx+gap:idx+gap+2,:]]), (4, 2))
            cells2Dsimplices.append(copy.deepcopy(tempSimp))
    
            for jdx in range(1, len(intervals)):
                tempSimp[:,1] += self.diameters[level-1]
                cells2Dsimplices.append(copy.deepcopy(tempSimp))    

        cells2Dsimplices = np.array(cells2Dsimplices)
        
        
        # Generate 3D mesh
        cubes = OrderedDict()
        cellCentroids = []
        cCount = 0
        for cell in cells2Dsimplices:
            
            for z in range(0, dims[2], self.diameters[level-1]):
                
                keyName = 'cell ' + str(cCount)
                singleCell = np.reshape(np.array([np.c_[cell, z*np.ones(len(cell))+offset[2]], np.c_[cell, (z+self.diameters[level-1])*np.ones(len(cell))+offset[2]]]), (8, 3))
                
                # Calculate centroids
                centroid = [np.sum(singleCell[:,0])/8.0, np.sum(singleCell[:,1])/8.0, np.sum(singleCell[:,2])/8.0]
                cellCentroids.append(copy.copy(centroid))
                
                cubes[keyName] = singleCell
                cCount += 1
        
        cellCoordinates = np.array(cubes.values())
        cubes.clear()
 
 
        return cellCentroids, cellCoordinates        
         
         
        
    def fill(self, level, node):
        '''Fill node and all branches of that node'''
        
        if level == 0: 
            return 'Root node cannot be filled.'
            
        self.tree[level][node][1] = True
        
        # Recursive fill        
        if self.tree[level][node][2]:
            return
        else:
            children = self.tree[level][node][4] 
            for child in children:
                self.tree[level+1][child][1] = True
                self.fill(level+1, child)
                
        
    def isFull(self, level, node):
        '''Is the node or any of its branches full
        [node number, full, leaf, parent, [children], coordinates, [centroid]]'''
        
        # Recursive search   
        if self.tree[level][node][1]: # if full
            return True
        else:
            if not self.tree[level][node][2]: # if not a leaf
                children = self.tree[level][node][4] 
                for child in children:
#                    print level+1, child
                    full = self.isFull(level+1, child)
                    if full:
                        return True
  
        
    def children(self, level, node):
        '''Node's children'''
        
        if level == 0:
            return self.tree[level][4]
            
        return self.tree[level][node][4]
        
        
    def parent(self, level, node):
        '''Node's parent'''
        
        if level == 0:
            return 'Root node has no parent.'
            
        return [level-1, self.tree[level][node][3]]
        

    def euclidean(self, x, y):
        '''Calculates Euclidean distance between 2 points'''
        
        x = np.array(x)
        y = np.array(y)
        return np.linalg.norm(x-y)


# Uncomment if using matplotlib & pyplot            
#    def plotCell(self, level, node, col='r', mark='o'):
#        '''Plot node corners'''
#        
#        ax.scatter(self.tree[level][node][5][:,0], self.tree[level][node][5][:,1], self.tree[level][node][5][:,2], c=col, marker=mark)
        
 
       
# Set up spheres in Abaqus, comment out if not running in Abaqus
def createAndPlaceSphere(diameter, centroid, sph):
    '''Creates & meshes spheres in the volume'''
    
    sphName = 'Sphere-' + str(sph)
    setName = 'SphereSet-' + str(sph)
    #rad = 0.1 + np.random.randint(10) * 0.1
    cent = diameter / 2.0
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=5.0)
    mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, -2.5), point2=(0.0, 2.5))
    mdb.models['Model-1'].sketches['__profile__'].FixedConstraint(entity=mdb.models['Model-1'].sketches['__profile__'].geometry[2])
    mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, -cent), point2=(0.0, cent)) #  direction=COUNTERCLOCKWISE,
    #print cent, rad
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name=sphName, type=DISCRETE_RIGID_SURFACE)
    mdb.models['Model-1'].parts[sphName].BaseShellRevolve(angle=360.0, flipRevolveDirection=OFF, sketch=mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']

    # Mesh spheres
    mdb.models['Model-1'].parts[sphName].setMeshControls(elemShape=QUAD, regions=mdb.models['Model-1'].parts[sphName].faces.getSequenceFromMask(('[#1 ]', ), ))
    mdb.models['Model-1'].parts[sphName].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=0.33)
    mdb.models['Model-1'].parts[sphName].generateMesh()
	
    # Create set
    mdb.models['Model-1'].parts[sphName].Set(faces=mdb.models['Model-1'].parts[sphName].faces.getSequenceFromMask(('[#1 ]', ), ), name=setName)
	
    # Reference points & inertia
    intName = 'Inertia-' + str(sph)
    refPointName = 'RP-' + str(sph)
    mdb.models['Model-1'].parts[sphName].ReferencePoint(point=(0.0, 0.0, 0.0))
    mdb.models['Model-1'].parts[sphName].Set(name=refPointName, referencePoints=(mdb.models['Model-1'].parts[sphName].referencePoints[5], ))
    volume = 0.75*np.pi*cent**3
    density = 2.33e-09
    mass = density / volume
    momInertia = (2.0/3.0)*mass*cent**3
    mdb.models['Model-1'].parts[sphName].engineeringFeatures.PointMassInertia(alpha=0.0, composite=0.0, i11=momInertia, i22=momInertia, i33=momInertia, mass=momInertia, name=intName, region=mdb.models['Model-1'].parts[sphName].sets[refPointName])
    
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=sphName, part=mdb.models['Model-1'].parts[sphName])
    xTrans = centroid[0]
    yTrans = centroid[1]
    zTrans = centroid[2]
#    mdb.models['Model-1'].rootAssembly.instances[sphName].translate(vector=(xTrans, yTrans, zTrans))
    mdb.models['Model-1'].rootAssembly.translate(instanceList=(sphName, ), vector=(xTrans, yTrans, zTrans))
    
def createAndCutSphere(diameter, centroid, sph):
    
    # Create sphere 
    sphName = 'Sphere-' + str(sph)          # Sphere thath is to be removed from root sphere
    radius = diameter/2.0
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=100.0)
    mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, -50.0), point2=(0.0, 50.0))
    mdb.models['Model-1'].sketches['__profile__'].FixedConstraint(entity=mdb.models['Model-1'].sketches['__profile__'].geometry[2])
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, radius), point2=(0.0, -radius))
    mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[3])
    mdb.models['Model-1'].sketches['__profile__'].ParallelConstraint(addUndoState=False, entity1=mdb.models['Model-1'].sketches['__profile__'].geometry[2], entity2=mdb.models['Model-1'].sketches['__profile__'].geometry[3])
    mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(addUndoState=False, entity1=mdb.models['Model-1'].sketches['__profile__'].vertices[0], entity2=mdb.models['Model-1'].sketches['__profile__'].geometry[2])
    mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(addUndoState=False, entity1=mdb.models['Model-1'].sketches['__profile__'].vertices[1], entity2=mdb.models['Model-1'].sketches['__profile__'].geometry[2])
    mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0), direction=CLOCKWISE, point1=(0.0, radius), point2=(0.0, -radius))
    mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(addUndoState=False, entity1=mdb.models['Model-1'].sketches['__profile__'].vertices[2], entity2=mdb.models['Model-1'].sketches['__profile__'].geometry[3])
    mdb.models['Model-1'].sketches['__profile__'].EqualDistanceConstraint(addUndoState=False, entity1=mdb.models['Model-1'].sketches['__profile__'].vertices[0], entity2=mdb.models['Model-1'].sketches['__profile__'].vertices[1], midpoint=mdb.models['Model-1'].sketches['__profile__'].vertices[2])
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name=sphName, type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts[sphName].BaseSolidRevolve(angle=360.0, flipRevolveDirection=OFF, sketch=mdb.models['Model-1'].sketches['__profile__'])
    
    # Place sphere to be cut instance    
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=sphName, part=mdb.models['Model-1'].parts[sphName])
    xTrans = centroid[0]
    yTrans = centroid[1]
    zTrans = centroid[2]
    mdb.models['Model-1'].rootAssembly.translate(instanceList=(sphName, ), vector=(xTrans, yTrans, zTrans))
    
    
def rootSphere(diam,dim, numberOfSpheres):
    
    rootName = 'rootSphere'
    radius = diam / 2.0
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=100.0)
    mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, -50.0), point2=(0.0, 50.0))
    mdb.models['Model-1'].sketches['__profile__'].FixedConstraint(entity=mdb.models['Model-1'].sketches['__profile__'].geometry[2])
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, radius), point2=(0.0, -radius))
    mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[3])
    mdb.models['Model-1'].sketches['__profile__'].ParallelConstraint(addUndoState=False, entity1=mdb.models['Model-1'].sketches['__profile__'].geometry[2], entity2=mdb.models['Model-1'].sketches['__profile__'].geometry[3])
    mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(addUndoState=False, entity1=mdb.models['Model-1'].sketches['__profile__'].vertices[0], entity2=mdb.models['Model-1'].sketches['__profile__'].geometry[2])
    mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(addUndoState=False, entity1=mdb.models['Model-1'].sketches['__profile__'].vertices[1], entity2=mdb.models['Model-1'].sketches['__profile__'].geometry[2])
    mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0), direction=CLOCKWISE, point1=(0.0, radius), point2=(0.0, -radius))
    mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(addUndoState=False, entity1=mdb.models['Model-1'].sketches['__profile__'].vertices[2], entity2=mdb.models['Model-1'].sketches['__profile__'].geometry[3])
    mdb.models['Model-1'].sketches['__profile__'].EqualDistanceConstraint(addUndoState=False, entity1=mdb.models['Model-1'].sketches['__profile__'].vertices[0], entity2=mdb.models['Model-1'].sketches['__profile__'].vertices[1], midpoint=mdb.models['Model-1'].sketches['__profile__'].vertices[2])
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name=sphName, type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts[sphName].BaseSolidRevolve(angle=360.0, flipRevolveDirection=OFF, sketch=mdb.models['Model-1'].sketches['__profile__'])
    
    # Add & cut instance    
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=sphName, part=mdb.models['Model-1'].parts[sphName])
    xTrans = dim [0]/ 2.0
    yTrans = dim [1]/ 2.0
    zTrans = dim [2]/ 2.0
    mdb.models['Model-1'].rootAssembly.translate(instanceList=(sphName, ), vector=(xTrans, yTrans, zTrans))
    
    # Cut operation - this will only work if script is created then run
    cutList = ''
    for sphNum in xrange(1, numberOfSpheres+1):
        cutList = cutList + 'mdb.models[\'Model-1\'].rootAssembly.instances[\'Sphere-' + str(sphNum) + '\'], '
    mdb.models['Model-1'].rootAssembly.InstanceFromBooleanCut(cuttingInstances=(mdb.models['Model-1'].rootAssembly.instances[sphName], ), instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances[rootSphere], name=newRoot, originalInstances=DELETE)

        

#####################################
#                                   #
#  Change these paramters as needed #
#                                   #   
#####################################

dimensions     = np.array([32, 32, 32]) # Global volume dimensions (does not have to be a cube)
diameters      = np.array([8, 2, 1])    # Increasing the number of diameters increases the set up time & memory usage
volumeFraction = 0.4                    # Target volume fraction, the higher this value the longer the execution, if too high it may never be attainable! 
isSphere       = False                  # If this flag is set to 'True' then fill space is treated as spherical & global dimensions need to represent a cube
                                        # If set to 'False' then the global dimensions give the shape of the volume
sponge         = False                  # Create negative space if set 'True'

myVol = VolumeFill(dimensions, diameters, isSphere) # This creates the tree structure that represents the volume

# Uncomment if using matplotlib & pyplot 
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.set_xlabel('X axis')
#ax.set_ylabel('Y axis')
#ax.set_zlabel('Z axis') 


# Populate space         

if isSphere:
    volume = (4.0/3.0) * np.pi * (np.max(dimensions)/2.0)**3
else:
    volume = dimensions[0]*dimensions[1]*dimensions[2]
    
fillCount = 0
sphVolume = 0
rootCentre = dimensions[0] / 2.0

# For FreeCAD
#doc = FreeCAD.newDocument("Sponge")
#rootSphere = Part.makeSphere(dimensions[0] / 2.0)
#rootTrans = Base.Vector(rootCentre, rootCentre, rootCentre)
#rootSphere.translate(rootTrans)

while (sphVolume/volume) < volumeFraction: 

    chosenLevel = np.random.randint(1, len(diameters)+1)
    chosenNode = np.random.randint(0, len(myVol.tree[chosenLevel]))

    
    if myVol.isFull(chosenLevel, chosenNode) == None:
        
        fillCount += 1
        myVol.fill(chosenLevel, chosenNode)
        
        if sponge:
            createAndCutSphere(diameters[chosenLevel-1], myVol.tree[chosenLevel][chosenNode][6], fillCount) # Un/comment for Abaqus
            sphere = Part.makeSphere((diameters[chosenLevel-1] + 0.5)/2.0)
            spongeCut = Base.Vector(myVol.tree[chosenLevel][chosenNode][6][0], myVol.tree[chosenLevel][chosenNode][6][1], myVol.tree[chosenLevel][chosenNode][6][2])
            sphere.translate(spongeCut)
            rootSphere = rootSphere.cut(sphere)
        else:
            createAndPlaceSphere(diameters[chosenLevel-1], myVol.tree[chosenLevel][chosenNode][6], fillCount) # Un/comment for Abaqus
            
        sphVolume += (4.0/3.0)*np.pi*(diameters[chosenLevel-1]/2.0)**3
        print('Volume fraction:', sphVolume/volume)
        
#if sponge:
#    rootSphere(dimensions[0], dimensions, fillCount)

#Part.show(rootSphere)

# __objs__=[]
# __objs__.append(FreeCAD.getDocument("Sponge").getObject("rootSphere"))
# Part.export(__objs__,u"C:/Users/User/Documents/FreeCAD scripts/Sponge_1.iges")

