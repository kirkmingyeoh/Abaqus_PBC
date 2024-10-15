# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:12:03 2024

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

Control point set names
V1 - Controls rigid body motion
V2 - Controls relative deformation along x-direction
V5 - Controls relative deformation along y-direction
V4 - Controls relative deformation along z-direction

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
ModelName = '3D_RVE_Test' # Name of CAE model
PartName = 'RVE' # Name of RVE Part
DP = 10 # For dealing with 1e-37 kind of numbers, only an issue when running through CAE

InstanceName = PartName + '-1'

### Define functions
def TakeVertexOut(face):
    del face[0]
    del face[-1]
    return face   

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

def SortListofNodes2D(faceN,coord1,coord2): # 0 for x; 1 for y; 2 for z; gives node numbers ready to be called with python (already -1)
    oldlistC1 = [] #first coordinate to be sorted along
    newlistN = []
    
    #Sort by first coodinate
    for i in range(len(faceN)):
        oldlistC1.append(RVENodalCoord[faceN[i]][coord1])
        
    newlistC1 = sorted(list(dict.fromkeys(oldlistC1)))
        
    for i in range(len(newlistC1)):
        C1 = newlistC1[i]
        sublistN = []
        sublistC2 = []
        for j in range(len(faceN)):
            C1N = RVENodalCoord[faceN[j]][coord1]
            
            if (C1N==C1):
                sublistN.append(faceN[j])
                sublistC2.append(RVENodalCoord[faceN[j]][coord2])
                
        newlistC2 = sorted(sublistC2)
        for j in range(len(sublistN)):
            Nindex = sublistC2.index(newlistC2[j])
            newlistN.append(sublistN[Nindex]) 
    
    return newlistN

def ExcludeNodes(faceN,coord1,coord2,coord3): # coord1, coord2 and coord3 are the coordinates to exclude, to be given as lists
    newlistN = []
    
    for i in range(len(faceN)):
        if RVENodalCoord[faceN[i]][0] not in coord1:
            if RVENodalCoord[faceN[i]][1] not in coord2:
                if RVENodalCoord[faceN[i]][2] not in coord3:
                    newlistN.append(faceN[i])
    
    return newlistN

### Obtaining and processing RVE information
# Obtaining all nodes of the RVE and their coordinates
RVENodalCoord = []
RVE_ListX = []
RVE_ListY = []
RVE_ListZ = []
Nnode = len(mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes)
for i in range(Nnode):
    x = round(mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[i].coordinates[0],DP)
    y = round(mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[i].coordinates[1],DP)
    z = round(mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[i].coordinates[2],DP)
    RVENodalCoord.append([x,y,z])
    RVE_ListX.append(x)
    RVE_ListY.append(y)
    RVE_ListZ.append(z)

# Finding the smallest and largest nodal coordinate in all 3 directions
xMin = min(RVE_ListX)
xMax = max(RVE_ListX)
yMin = min(RVE_ListY)
yMax = max(RVE_ListY)
zMin = min(RVE_ListZ)
zMax = max(RVE_ListZ)
del RVE_ListX
del RVE_ListY
del RVE_ListZ
    
# Sorting the RVE nodes
FaceLNodes,FaceRNodes,FaceBaNodes,FaceFNodes,FaceBNodes,FaceTNodes = [],[],[],[],[],[] # Faces along x, y and z, 2 each
EdgeLNodes,EdgeRNodes,EdgeTNodes,EdgeBNodes = [],[],[],[] # Edges along FaceBa
Edge1Nodes,Edge2Nodes,Edge3Nodes,Edge4Nodes = [],[],[],[] # Edges parallel to z direction, between FaceBa and FaceF
for i in range(len(RVENodalCoord)):
    if RVENodalCoord[i][0] == xMin:
        FaceLNodes.append(i)
    if RVENodalCoord[i][0] == xMax:
        FaceRNodes.append(i)
    if RVENodalCoord[i][1] == yMin:
        FaceBaNodes.append(i)
    if RVENodalCoord[i][1] == yMax:
        FaceFNodes.append(i)
    if RVENodalCoord[i][2] == zMin:
        FaceBNodes.append(i)
    if RVENodalCoord[i][2] == zMax:
        FaceTNodes.append(i)
for i in range(len(FaceBaNodes)):
    if FaceBaNodes[i] in FaceLNodes:
        EdgeLNodes.append(FaceBaNodes[i])
        if FaceBaNodes[i] in FaceBNodes:
            V1 = FaceBaNodes[i]
        if FaceBaNodes[i] in FaceTNodes:
            V4 = FaceBaNodes[i]
    if FaceBaNodes[i] in FaceRNodes:
        EdgeRNodes.append(FaceBaNodes[i])
        if FaceBaNodes[i] in FaceBNodes:
            V2 = FaceBaNodes[i]
        if FaceBaNodes[i] in FaceTNodes:
            V3 = FaceBaNodes[i]
    if FaceBaNodes[i] in FaceBNodes:
        EdgeBNodes.append(FaceBaNodes[i])
    if FaceBaNodes[i] in FaceTNodes:
        EdgeTNodes.append(FaceBaNodes[i])
for i in range(len(FaceBNodes)):
    if FaceBNodes[i] in FaceLNodes:
        Edge1Nodes.append(FaceBNodes[i])
        if FaceBNodes[i] in FaceFNodes:
            V5 = FaceBNodes[i]
    if FaceBNodes[i] in FaceRNodes:
        Edge2Nodes.append(FaceBNodes[i])
        if FaceBNodes[i] in FaceFNodes:
            V6 = FaceBNodes[i]
for i in range(len(FaceTNodes)):
    if FaceTNodes[i] in FaceLNodes:
        Edge4Nodes.append(FaceTNodes[i])
        if FaceTNodes[i] in FaceFNodes:
            V8 = FaceTNodes[i]
    if FaceTNodes[i] in FaceRNodes:
        Edge3Nodes.append(FaceTNodes[i])
        if FaceTNodes[i] in FaceFNodes:
            V7 = FaceTNodes[i]

### Pairing nodes
EdgeLNodes = TakeVertexOut(SortListofNodes1D(EdgeLNodes,2))
EdgeRNodes = TakeVertexOut(SortListofNodes1D(EdgeRNodes,2)) 
EdgeBNodes = TakeVertexOut(SortListofNodes1D(EdgeBNodes,0))
EdgeTNodes = TakeVertexOut(SortListofNodes1D(EdgeTNodes,0))
    
Edge1Nodes = TakeVertexOut(SortListofNodes1D(Edge1Nodes,1))
Edge2Nodes = TakeVertexOut(SortListofNodes1D(Edge2Nodes,1))
Edge3Nodes = TakeVertexOut(SortListofNodes1D(Edge3Nodes,1))
Edge4Nodes = TakeVertexOut(SortListofNodes1D(Edge4Nodes,1))

# Sort and pair all nodes on FaceBa and FaceF
FaceBaNodes = SortListofNodes2D(FaceBaNodes,0,2)
FaceFNodes = SortListofNodes2D(FaceFNodes,0,2)
PairingFacesBaF = []
for i in range(len(FaceBaNodes)):
    # Skip V1 and V5 which are control nodes
    if FaceBaNodes[i] == V1:
        continue
    
    Temp = []
    Temp.append(FaceBaNodes[i])
    x1 = RVENodalCoord[FaceBaNodes[i]][0]
    z1 = RVENodalCoord[FaceBaNodes[i]][2]
    
    Dist = []
    for j in range(len(FaceFNodes)):
        x2 = RVENodalCoord[FaceFNodes[j]][0]
        z2 = RVENodalCoord[FaceFNodes[j]][2]
        
        Dist.append(math.sqrt(pow(x2-x1,2)+pow(z2-z1,2)))
        
    N = Dist.index(min(Dist))
    Temp.append(FaceFNodes[N])
    FaceFNodes.pop(N)
    PairingFacesBaF.append(Temp)    

# Sort and pair all nodes on FaceL and FaceR not on their edges and vertices
# Exclude out nodes on the edges and vertices before sorting
FaceLNodes = SortListofNodes2D(ExcludeNodes(FaceLNodes,[],[RVENodalCoord[V3][1],RVENodalCoord[V6][1]],[RVENodalCoord[V3][2],RVENodalCoord[V6][2]]),1,2)
FaceRNodes = SortListofNodes2D(ExcludeNodes(FaceRNodes,[],[RVENodalCoord[V3][1],RVENodalCoord[V6][1]],[RVENodalCoord[V3][2],RVENodalCoord[V6][2]]),1,2)
PairingFacesLR = []
for i in range(len(FaceLNodes)):
    Temp = []
    Temp.append(FaceLNodes[i])
    y1 = RVENodalCoord[FaceLNodes[i]][1]
    z1 = RVENodalCoord[FaceLNodes[i]][2]
    
    Dist = []
    for j in range(len(FaceRNodes)):
        y2 = RVENodalCoord[FaceRNodes[j]][1]
        z2 = RVENodalCoord[FaceRNodes[j]][2]
        
        Dist.append(math.sqrt(pow(y2-y1,2)+pow(z2-z1,2)))
        
    N = Dist.index(min(Dist))
    Temp.append(FaceRNodes[N])
    FaceRNodes.pop(N)
    PairingFacesLR.append(Temp)

# Sort and pair all nodes on FaceB and FaceT not on their edges and vertices
# Exclude out nodes on the edges and vertices before sorting
FaceBNodes = SortListofNodes2D(ExcludeNodes(FaceBNodes,[RVENodalCoord[V3][0],RVENodalCoord[V8][0]],[RVENodalCoord[V3][1],RVENodalCoord[V8][1]],[]),0,1)
FaceTNodes = SortListofNodes2D(ExcludeNodes(FaceTNodes,[RVENodalCoord[V3][0],RVENodalCoord[V8][0]],[RVENodalCoord[V3][1],RVENodalCoord[V8][1]],[]),0,1)
PairingFacesBT = []
for i in range(len(FaceBNodes)):
    Temp = []
    Temp.append(FaceBNodes[i])
    x1 = RVENodalCoord[FaceBNodes[i]][0]
    y1 = RVENodalCoord[FaceBNodes[i]][1]
    
    Dist = []
    for j in range(len(FaceTNodes)):
        x2 = RVENodalCoord[FaceTNodes[j]][0]
        y2 = RVENodalCoord[FaceTNodes[j]][1]
        
        Dist.append(math.sqrt(pow(x2-x1,2)+pow(y2-y1,2)))
        
    N = Dist.index(min(Dist))
    Temp.append(FaceTNodes[N])
    FaceTNodes.pop(N)
    PairingFacesBT.append(Temp)

### Calling sets and setting up PBCs
# Calling sets for corner nodes        
mdb.models[ModelName].rootAssembly.Set(name='V1',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V1:V1+1])
mdb.models[ModelName].rootAssembly.Set(name='V2',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V2:V2+1])
mdb.models[ModelName].rootAssembly.Set(name='V3',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V3:V3+1])
mdb.models[ModelName].rootAssembly.Set(name='V4',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V4:V4+1])
mdb.models[ModelName].rootAssembly.Set(name='V5',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V5:V5+1])
mdb.models[ModelName].rootAssembly.Set(name='V6',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V6:V6+1])
mdb.models[ModelName].rootAssembly.Set(name='V7',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V7:V7+1])
mdb.models[ModelName].rootAssembly.Set(name='V8',nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[V8:V8+1])

# EdgeL and EdgeR
for i in range(len(EdgeLNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-EdgeLNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[EdgeLNodes[i]:EdgeLNodes[i]+1])
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-EdgeRNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[EdgeRNodes[i]:EdgeRNodes[i]+1])
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-EdgeLR-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-EdgeRNode'+str(i+1),j+1),
            (1.0,str(InstanceName)+'-EdgeLNode'+str(i+1),j+1),
            (1.0,'V2',j+1),
            (-1.0,'V1',j+1)))

# EdgeB and EdgeT
for i in range(len(EdgeBNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-EdgeBNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[EdgeBNodes[i]:EdgeBNodes[i]+1])
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-EdgeTNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[EdgeTNodes[i]:EdgeTNodes[i]+1])
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-EdgeBT-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-EdgeTNode'+str(i+1),j+1),
            (1.0,str(InstanceName)+'-EdgeBNode'+str(i+1),j+1),
            (1.0,'V4',j+1),
            (-1.0,'V1',j+1)))

# V3
for i in range(3):
    mdb.models[ModelName].Equation(name=str(InstanceName)+'-V3-DOF-'+str(i+1),terms=(
        (-1.0,'V3',i+1),
        (1.0,'V2',i+1),
        (1.0,'V4',i+1),
        (-1.0,'V1',i+1)))

# Edge1, Edge2, Edge3, Edge4
for i in range(len(Edge1Nodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-Edge1Node'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[Edge1Nodes[i]:Edge1Nodes[i]+1])
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-Edge2Node'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[Edge2Nodes[i]:Edge2Nodes[i]+1])
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-Edge3Node'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[Edge3Nodes[i]:Edge3Nodes[i]+1])
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-Edge4Node'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[Edge4Nodes[i]:Edge4Nodes[i]+1])
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-Edge14-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-Edge4Node'+str(i+1),j+1),
            (1.0,str(InstanceName)+'-Edge1Node'+str(i+1),j+1),
            (1.0,'V4',j+1),
            (-1.0,'V1',j+1)))
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-Edge23-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-Edge3Node'+str(i+1),j+1),
            (1.0,str(InstanceName)+'-Edge2Node'+str(i+1),j+1),
            (1.0,'V4',j+1),
            (-1.0,'V1',j+1)))
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-Edge12-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-Edge2Node'+str(i+1),j+1),
            (1.0,str(InstanceName)+'-Edge1Node'+str(i+1),j+1),
            (1.0,'V2',j+1),
            (-1.0,'V1',j+1)))

# FaceL and FaceR
for i in range(len(PairingFacesLR)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceLNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[PairingFacesLR[i][0]:PairingFacesLR[i][0]+1])
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceRNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[PairingFacesLR[i][1]:PairingFacesLR[i][1]+1])
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-FaceLR-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-FaceRNode'+str(i+1),j+1),
            (1.0,str(InstanceName)+'-FaceLNode'+str(i+1),j+1),
            (1.0,'V2',j+1),
            (-1.0,'V1',j+1)))

# FaceB and FaceT
for i in range(len(PairingFacesBT)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceBNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[PairingFacesBT[i][0]:PairingFacesBT[i][0]+1])
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceTNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[PairingFacesBT[i][1]:PairingFacesBT[i][1]+1])
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-FaceBT-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-FaceTNode'+str(i+1),j+1),
            (1.0,str(InstanceName)+'-FaceBNode'+str(i+1),j+1),
            (1.0,'V4',j+1),
            (-1.0,'V1',j+1)))

# FaceBa and FaceF
for i in range(len(PairingFacesBaF)):    
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceBaNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[PairingFacesBaF[i][0]:PairingFacesBaF[i][0]+1])
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceFNode'+str(i+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[PairingFacesBaF[i][1]:PairingFacesBaF[i][1]+1])
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-FaceBaF-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-FaceFNode'+str(i+1),j+1),
            (1.0,str(InstanceName)+'-FaceBaNode'+str(i+1),j+1),
            (1.0,'V5',j+1),
            (-1.0,'V1',j+1)))

