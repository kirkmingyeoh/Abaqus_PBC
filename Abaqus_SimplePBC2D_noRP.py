# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 14:21:52 2022

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

Control point set names
V1 - Controls rigid body motion
V2 - Controls relative deformation along x-direction
V4 - Controls relative deformation along y-direction

Apply the BCs using the respective set names

"""

import string
import math
import sys
import os
import numpy as np

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

### Update the following parameters
ModelName = '2D_RVE_Test' # Name of CAE model
PartName = 'RVE' # Name of RVE Part
DP = 10 # For dealing with 1e-37 kind of numbers, only an issue when running through CAE

InstanceName = PartName + '-1'

### Define functions
def SortListofNodes1D(faceN,coordinate): # 0 for x; 1 for y; 2 for z; gives node numbers ready to be called with python (already -1)
    newlist = []
    oldlist = []
    for i in range(len(faceN)):
        oldlist.append(RVENodalCoord[faceN[i]][coordinate])
    
    orderedlist = sorted(oldlist)
    for j in range(len(orderedlist)):
        ind = oldlist.index(orderedlist[j])
        newlist.append(faceN[ind])
    
    return newlist

def TakeVertexOut(face):
    del face[0]
    del face[-1]
    return face 

### Obtaining and processing RVE information
# Obtaining all nodes of the RVE and their coordinates
RVENodalCoord = []
RVE_ListX = []
RVE_ListY = []
Nnode = len(mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes)
for i in range(Nnode):
    x = round(mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[i].coordinates[0],DP)
    y = round(mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[i].coordinates[1],DP)
    RVENodalCoord.append([x,y])
    RVE_ListX.append(x)
    RVE_ListY.append(y)

# Finding the smallest and largest nodal coordinate in all 3 directions
xMin = min(RVE_ListX)
xMax = max(RVE_ListX)
yMin = min(RVE_ListY)
yMax = max(RVE_ListY)
del RVE_ListX
del RVE_ListY

# Sorting the RVE nodes
FaceLNodes,FaceRNodes,FaceBNodes,FaceTNodes = [],[],[],[]
for i in range(len(RVENodalCoord)):
    if RVENodalCoord[i][0] == xMin:
        FaceLNodes.append(i)
    if RVENodalCoord[i][0] == xMax:
        FaceRNodes.append(i)
    if RVENodalCoord[i][1] == yMin:
        FaceBNodes.append(i)
    if RVENodalCoord[i][1] == yMax:
        FaceTNodes.append(i)
for i in range(len(FaceLNodes)):
    if FaceLNodes[i] in FaceBNodes:
        V1 = FaceLNodes[i]
    if FaceLNodes[i] in FaceTNodes:
        V4 = FaceLNodes[i]
for i in range(len(FaceRNodes)):
    if FaceRNodes[i] in FaceBNodes:
        V2 = FaceRNodes[i]
    if FaceRNodes[i] in FaceTNodes:
        V3 = FaceRNodes[i]

# Sort all nodes on boundary faces but not on the vertices
FaceLNodes = TakeVertexOut(SortListofNodes1D(FaceLNodes,1))
FaceRNodes = TakeVertexOut(SortListofNodes1D(FaceRNodes,1))
FaceBNodes = TakeVertexOut(SortListofNodes1D(FaceBNodes,0))
FaceTNodes = TakeVertexOut(SortListofNodes1D(FaceTNodes,0))

### Calling sets and setting up PBCs
# Calling sets for corner nodes 
mdb.models[ModelName].rootAssembly.Set(name='V1',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V1:V1+1])
mdb.models[ModelName].rootAssembly.Set(name='V2',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V2:V2+1])
mdb.models[ModelName].rootAssembly.Set(name='V3',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V3:V3+1])
mdb.models[ModelName].rootAssembly.Set(name='V4',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V4:V4+1])

# FaceL and FaceR
for i in range(len(FaceLNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceLNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[FaceLNodes[i]:FaceLNodes[i]+1])
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceRNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[FaceRNodes[i]:FaceRNodes[i]+1])
    
    for j in range(2):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-FaceLR-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-FaceRNode'+str(i+1),j+1),
            (1.0,str(InstanceName)+'-FaceLNode'+str(i+1),j+1),
            (1.0,'V2',j+1),
            (-1.0,'V1',j+1)))
    
# FaceB and FaceT
for i in range(len(FaceBNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceBNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[FaceBNodes[i]:FaceBNodes[i]+1])
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceTNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[FaceTNodes[i]:FaceTNodes[i]+1])
    
    for j in range(2):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-FaceBT-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-FaceTNode'+str(i+1),j+1),
            (1.0,str(InstanceName)+'-FaceBNode'+str(i+1),j+1),
            (1.0,'V4',j+1),
            (-1.0,'V1',j+1)))

# V3
for i in range(2):
    mdb.models[ModelName].Equation(name=str(InstanceName)+'-V3-DOF-'+str(i+1),terms=(
        (-1.0,'V3',i+1),
        (1.0,'V2',i+1),
        (1.0,'V4',i+1),
        (-1.0,'V1',i+1)))












