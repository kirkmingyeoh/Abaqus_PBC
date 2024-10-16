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

**Note: This script assumes CPS4 or CPS4R elements are used to mesh the RVE

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
def TakeVertexOut(edge):
    del edge[0]
    del edge[-1]
    return edge   

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

def SortFaceConnect(Face1Ele,Face1Nodes,C1,C2): # Sorts connectivity of nodes of elements on the master faces of the RVE (FaceL, FaceBa and FaceB)
    Face1Connect = [[] for x in range(len(Face1Ele))]
    
    for i in range(len(Face1Ele)):
        TempConnect = []
        for j in range(8):
            Node = RVENodalConnect[Face1Ele[i]][j]
            if (Node in Face1Nodes):
                TempConnect.append(Node)               
        Y = []
        Z = []
        Ang = []
        for j in range(4):
            Y.append(RVENodalCoord[TempConnect[j]][C1])
            Z.append(RVENodalCoord[TempConnect[j]][C2])          
        y0 = sum(Y)/4
        z0 = sum(Z)/4        
        for j in range(4):
            y1 = Y[j]-y0
            z1 = Z[j]-z0
            if z1 == 0:
                if y1 > 0:
                    theta = math.pi*0.5
                else:
                    theta = math.pi*1.5
            else:
                theta = math.atan(y1/z1)
            if (z1 < 0):
                theta = theta+math.pi            
            if (theta < 0):
                theta = theta+2*(math.pi)           
            Ang.append(theta*360/(2*(math.pi)))            
        SAng = sorted(Ang)
        
        for j in range(4):
            N = Ang.index(SAng[j])
            Face1Connect[i].append(TempConnect[N])           
        Face1Connect[i][0],Face1Connect[i][1],Face1Connect[i][2],Face1Connect[i][3] = Face1Connect[i][2],Face1Connect[i][3],Face1Connect[i][0],Face1Connect[i][1]
        
    return Face1Connect

def FaceNodetoElePairing(Face2Nodes,Face1Nodes,Face1Connect,C1,C2): #2 is S, 1 is M for nodes and elements; C1 and C2 are the 2 in-plane directions 
    Face1Pairing = [[] for x in range(len(Face2Nodes))]
    for i in range(len(Face2Nodes)):
        yR = RVENodalCoord[Face2Nodes[i]][C1]
        zR = RVENodalCoord[Face2Nodes[i]][C2]
        
        Dist = []
        for j in range(len(Face1Nodes)):
            yL = RVENodalCoord[Face1Nodes[j]][C1]
            zL = RVENodalCoord[Face1Nodes[j]][C2]
            Dist.append(math.sqrt(pow(yR-yL,2)+pow(zR-zL,2)))
            
        Node1 = Face1Nodes[Dist.index(min(Dist))]
        
        Ele1 = []
        for j in range(len(Face1Connect)):
            if Node1 in Face1Connect[j]:
                Ele1.append(Face1Connect[j])
                
        for j in range(len(Ele1)):
            y1 = RVENodalCoord[Ele1[j][0]][C1]
            y2 = RVENodalCoord[Ele1[j][1]][C1]
            y3 = RVENodalCoord[Ele1[j][2]][C1]
            y4 = RVENodalCoord[Ele1[j][3]][C1]
            z1 = RVENodalCoord[Ele1[j][0]][C2]
            z2 = RVENodalCoord[Ele1[j][1]][C2]
            z3 = RVENodalCoord[Ele1[j][2]][C2]
            z4 = RVENodalCoord[Ele1[j][3]][C2]
            
            C = np.array([[1,y1,z1,y1*z1],[1,y2,z2,y2*z2],[1,y3,z3,y3*z3],[1,y4,z4,y4*z4]])
            [a0,a1,a2,a3] = np.transpose(np.dot(np.linalg.inv(C),np.transpose(np.array([-1,1,1,-1]))))
            [b0,b1,b2,b3] = np.transpose(np.dot(np.linalg.inv(C),np.transpose(np.array([-1,-1,1,1]))))
            
            tsi = round((a0 + a1*yR + a2*zR + a3*yR*zR),DP)
            eta = round((b0 + b1*yR + b2*zR + b3*yR*zR),DP)
            
            if (tsi<=1) and (tsi>=-1) and (eta<=1) and (eta>=-1):
                Face1Pairing[i].append(Ele1[j])
                SF = [0.25*(1-eta)*(1-tsi),0.25*(1-eta)*(1+tsi),0.25*(1+eta)*(1+tsi),0.25*(1+eta)*(1-tsi)]
                Face1Pairing[i].append(SF)
                break
        
    return Face1Pairing

def EdgeNodetoNodePairing(Edge2Nodes,Edge1Nodes,C1): #2 is S, 1 is M; C1 is the direction along the edge
    PairingEdgeNodes = []
    
    Start = 0
    for i in range(len(Edge2Nodes)):
        Coord1 = RVENodalCoord[Edge2Nodes[i]][C1]
        for j in range(Start,len(Edge1Nodes)):
            Coord2 = RVENodalCoord[Edge1Nodes[j]][C1]
            if Coord2 == Coord1:
                PairingEdgeNodes.append([[Edge1Nodes[j],1.0]])
                Start = j
                break
            if Coord2 > Coord1:
                Diff = Coord2 - Coord1
                Frac = Diff/(Coord2 - RVENodalCoord[Edge1Nodes[j-1]][C1])
                PairingEdgeNodes.append([[Edge1Nodes[j-1],Frac],[Edge1Nodes[j],1.0-Frac]])
                break
    
    return PairingEdgeNodes #[node sorted position in list, 1.0] or [[node1 sorted position in list, Frac],[node2 sorted position in list, Frac]]

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
        
# Obtaining all elements of the RVE and their connectivities
RVENodalConnect = []
Nele = len(mdb.models[ModelName].rootAssembly.instances[InstanceName].elements)
for i in range(Nele):
    Temp = []
    for j in range(8):
        Temp.append(mdb.models[ModelName].rootAssembly.instances[InstanceName].elements[i].connectivity[j])
    RVENodalConnect.append(Temp)

# Finding elements located on the master faces (FaceL, FaceBa, FaceB)
FaceLEle, FaceBaEle, FaceBEle = [],[],[]
for i in range(len(RVENodalConnect)):
    for j in range(8):
        N = RVENodalConnect[i][j]
        if N in FaceLNodes:
            FaceLEle.append(i)
        if N in FaceBaNodes:
            FaceBaEle.append(i)
        if N in FaceBNodes:
            FaceBEle.append(i)
        
### Pairing nodes
# Sort and pair nodes on FaceR to node group on FaceL
# Exclude out nodes on the edges and vertices of FaceR before sorting 
FaceRNodes = SortListofNodes2D(ExcludeNodes(FaceRNodes,[],[RVENodalCoord[V3][1],RVENodalCoord[V6][1]],[RVENodalCoord[V3][2],RVENodalCoord[V6][2]]),1,2)
FaceLConnect = SortFaceConnect(FaceLEle,FaceLNodes,1,2)
PairingFacesLR = FaceNodetoElePairing(FaceRNodes,FaceLNodes,FaceLConnect,1,2) 

# Sort and pair nodes on FaceT to node group on FaceB
# Exclude out nodes on the edges and vertices of FaceT before sorting 
FaceTNodes = SortListofNodes2D(ExcludeNodes(FaceTNodes,[RVENodalCoord[V3][0],RVENodalCoord[V8][0]],[RVENodalCoord[V3][1],RVENodalCoord[V8][1]],[]),0,1)
FaceBConnect = SortFaceConnect(FaceBEle,FaceBNodes,0,1)
PairingFacesBT = FaceNodetoElePairing(FaceTNodes,FaceBNodes,FaceBConnect,0,1) 

# Sort and pair nodes on FaceF to node group on FaceBa
FaceFNodes = SortListofNodes2D(FaceFNodes,0,2)
FaceBaConnect = SortFaceConnect(FaceBaEle,FaceBaNodes,0,2)
PairingFacesBaF = FaceNodetoElePairing(FaceFNodes,FaceBaNodes,FaceBaConnect,0,2) 

# Sorting the edge nodes, preserving corners on the master edges (EdgeL, EdgeB and Edge1)
EdgeLNodes = SortListofNodes1D(EdgeLNodes,2)
EdgeBNodes = SortListofNodes1D(EdgeBNodes,0)
Edge1Nodes = SortListofNodes1D(Edge1Nodes,1)

EdgeRNodes = TakeVertexOut(SortListofNodes1D(EdgeRNodes,2))
EdgeTNodes = TakeVertexOut(SortListofNodes1D(EdgeTNodes,0))
Edge2Nodes = TakeVertexOut(SortListofNodes1D(Edge2Nodes,1))
Edge3Nodes = TakeVertexOut(SortListofNodes1D(Edge3Nodes,1))
Edge4Nodes = TakeVertexOut(SortListofNodes1D(Edge4Nodes,1))

# Pairing the nodes on the edges
PairingEdgesLR = EdgeNodetoNodePairing(EdgeRNodes,EdgeLNodes,2)
PairingEdgesBT = EdgeNodetoNodePairing(EdgeTNodes,EdgeBNodes,0)
PairingEdges14 = EdgeNodetoNodePairing(Edge4Nodes,Edge1Nodes,1)
PairingEdges13 = EdgeNodetoNodePairing(Edge3Nodes,Edge1Nodes,1)
PairingEdges12 = EdgeNodetoNodePairing(Edge2Nodes,Edge1Nodes,1)

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

# V3
for i in range(3):
    mdb.models[ModelName].Equation(name=str(InstanceName)+'-V3-DOF-'+str(i+1),terms=(
        (-1.0,'V3',i+1),
        (1.0,'V2',i+1),
        (1.0,'V4',i+1),
        (-1.0,'V1',i+1)))
    
# EdgeL and EdgeR
for i in range(len(EdgeLNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-EdgeLNode'+str(EdgeLNodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[EdgeLNodes[i]:EdgeLNodes[i]+1])
    
for i in range(len(EdgeRNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-EdgeRNode'+str(EdgeRNodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[EdgeRNodes[i]:EdgeRNodes[i]+1])
    
    if len(PairingEdgesLR[i]) == 1:
        for j in range(3):
            mdb.models[ModelName].Equation(name=str(InstanceName)+'-EdgeLR-'+str(i+1)+'-DOF-'+str(j+1),terms=(
                (-1.0,str(InstanceName)+'-EdgeRNode'+str(EdgeRNodes[i]+1),j+1),
                (1.0,str(InstanceName)+'-EdgeLNode'+str(PairingEdgesLR[i][0][0]+1),j+1),
                (1.0,'V2',j+1),
                (-1.0,'V1',j+1)))
    else:
        for j in range(3):
            mdb.models[ModelName].Equation(name=str(InstanceName)+'-EdgeLR-'+str(i+1)+'-DOF-'+str(j+1),terms=(
                (-1.0,str(InstanceName)+'-EdgeRNode'+str(EdgeRNodes[i]+1),j+1),
                (PairingEdgesLR[i][0][1],str(InstanceName)+'-EdgeLNode'+str(PairingEdgesLR[i][0][0]+1),j+1),
                (PairingEdgesLR[i][1][1],str(InstanceName)+'-EdgeLNode'+str(PairingEdgesLR[i][1][0]+1),j+1),
                (1.0,'V2',j+1),
                (-1.0,'V1',j+1)))

# EdgeB and EdgeT
for i in range(len(EdgeBNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-EdgeBNode'+str(EdgeBNodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[EdgeBNodes[i]:EdgeBNodes[i]+1])
    
for i in range(len(EdgeTNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-EdgeTNode'+str(EdgeTNodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[EdgeTNodes[i]:EdgeTNodes[i]+1])
    
    if len(PairingEdgesBT[i]) == 1:
        for j in range(3):
            mdb.models[ModelName].Equation(name=str(InstanceName)+'-EdgeBT-'+str(i+1)+'-DOF-'+str(j+1),terms=(
                (-1.0,str(InstanceName)+'-EdgeTNode'+str(EdgeTNodes[i]+1),j+1),
                (1.0,str(InstanceName)+'-EdgeBNode'+str(PairingEdgesBT[i][0][0]+1),j+1),
                (1.0,'V4',j+1),
                (-1.0,'V1',j+1)))
    else:
        for j in range(3):
            mdb.models[ModelName].Equation(name=str(InstanceName)+'-EdgeBT-'+str(i+1)+'-DOF-'+str(j+1),terms=(
                (-1.0,str(InstanceName)+'-EdgeTNode'+str(EdgeTNodes[i]+1),j+1),
                (PairingEdgesBT[i][0][1],str(InstanceName)+'-EdgeBNode'+str(PairingEdgesBT[i][0][0]+1),j+1),
                (PairingEdgesBT[i][1][1],str(InstanceName)+'-EdgeBNode'+str(PairingEdgesBT[i][1][0]+1),j+1),
                (1.0,'V4',j+1),
                (-1.0,'V1',j+1)))

# Edge1, Edge2, Edge3 and Edge4
for i in range(len(Edge1Nodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-Edge1Node'+str(Edge1Nodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[Edge1Nodes[i]:Edge1Nodes[i]+1])
    
for i in range(len(Edge4Nodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-Edge4Node'+str(Edge4Nodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[Edge4Nodes[i]:Edge4Nodes[i]+1])
    
    if len(PairingEdges14[i]) == 1:
        for j in range(3):
            mdb.models[ModelName].Equation(name=str(InstanceName)+'-Edge14-'+str(i+1)+'-DOF-'+str(j+1),terms=(
                (-1.0,str(InstanceName)+'-Edge4Node'+str(Edge4Nodes[i]+1),j+1),
                (1.0,str(InstanceName)+'-Edge1Node'+str(PairingEdges14[i][0][0]+1),j+1),
                (1.0,'V4',j+1),
                (-1.0,'V1',j+1)))
    else:
        for j in range(3):
            mdb.models[ModelName].Equation(name=str(InstanceName)+'-Edge14-'+str(i+1)+'-DOF-'+str(j+1),terms=(
                (-1.0,str(InstanceName)+'-Edge4Node'+str(Edge4Nodes[i]+1),j+1),
                (PairingEdges14[i][0][1],str(InstanceName)+'-Edge1Node'+str(PairingEdges14[i][0][0]+1),j+1),
                (PairingEdges14[i][1][1],str(InstanceName)+'-Edge1Node'+str(PairingEdges14[i][1][0]+1),j+1),
                (1.0,'V4',j+1),
                (-1.0,'V1',j+1)))

for i in range(len(Edge3Nodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-Edge3Node'+str(Edge3Nodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[Edge3Nodes[i]:Edge3Nodes[i]+1])
    
    if len(PairingEdges13[i]) == 1:
        for j in range(3):
            mdb.models[ModelName].Equation(name=str(InstanceName)+'-Edge13-'+str(i+1)+'-DOF-'+str(j+1),terms=(
                (-1.0,str(InstanceName)+'-Edge3Node'+str(Edge3Nodes[i]+1),j+1),
                (1.0,str(InstanceName)+'-Edge1Node'+str(PairingEdges13[i][0][0]+1),j+1),
                (1.0,'V3',j+1),
                (-1.0,'V1',j+1)))
    else:
        for j in range(3):
            mdb.models[ModelName].Equation(name=str(InstanceName)+'-Edge13-'+str(i+1)+'-DOF-'+str(j+1),terms=(
                (-1.0,str(InstanceName)+'-Edge3Node'+str(Edge3Nodes[i]+1),j+1),
                (PairingEdges13[i][0][1],str(InstanceName)+'-Edge1Node'+str(PairingEdges13[i][0][0]+1),j+1),
                (PairingEdges13[i][1][1],str(InstanceName)+'-Edge1Node'+str(PairingEdges13[i][1][0]+1),j+1),
                (1.0,'V3',j+1),
                (-1.0,'V1',j+1)))

for i in range(len(Edge2Nodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-Edge2Node'+str(Edge2Nodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[Edge2Nodes[i]:Edge2Nodes[i]+1])
    
    if len(PairingEdges12[i]) == 1:
        for j in range(3):
            mdb.models[ModelName].Equation(name=str(InstanceName)+'-Edge12-'+str(i+1)+'-DOF-'+str(j+1),terms=(
                (-1.0,str(InstanceName)+'-Edge2Node'+str(Edge2Nodes[i]+1),j+1),
                (1.0,str(InstanceName)+'-Edge1Node'+str(PairingEdges12[i][0][0]+1),j+1),
                (1.0,'V2',j+1),
                (-1.0,'V1',j+1)))
    else:
        for j in range(3):
            mdb.models[ModelName].Equation(name=str(InstanceName)+'-Edge12-'+str(i+1)+'-DOF-'+str(j+1),terms=(
                (-1.0,str(InstanceName)+'-Edge2Node'+str(Edge2Nodes[i]+1),j+1),
                (PairingEdges12[i][0][1],str(InstanceName)+'-Edge1Node'+str(PairingEdges12[i][0][0]+1),j+1),
                (PairingEdges12[i][1][1],str(InstanceName)+'-Edge1Node'+str(PairingEdges12[i][1][0]+1),j+1),
                (1.0,'V2',j+1),
                (-1.0,'V1',j+1)))

# FaceL and FaceR
for i in range(len(FaceLNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceLNode'+str(FaceLNodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[FaceLNodes[i]:FaceLNodes[i]+1])
    
for i in range(len(FaceRNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceRNode'+str(FaceRNodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[FaceRNodes[i]:FaceRNodes[i]+1])
    
    LNodes = PairingFacesLR[i][0]
    LSF = PairingFacesLR[i][1]
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-FaceLR-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-FaceRNode'+str(FaceRNodes[i]+1),j+1),
            (LSF[0],str(InstanceName)+'-FaceLNode'+str(LNodes[0]+1),j+1),
            (LSF[1],str(InstanceName)+'-FaceLNode'+str(LNodes[1]+1),j+1),
            (LSF[2],str(InstanceName)+'-FaceLNode'+str(LNodes[2]+1),j+1),
            (LSF[3],str(InstanceName)+'-FaceLNode'+str(LNodes[3]+1),j+1),
            (1.0,'V2',j+1),
            (-1.0,'V1',j+1)))

# FaceB and FaceT
for i in range(len(FaceBNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceBNode'+str(FaceBNodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[FaceBNodes[i]:FaceBNodes[i]+1])
    
for i in range(len(FaceTNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceTNode'+str(FaceTNodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[FaceTNodes[i]:FaceTNodes[i]+1])
    
    BNodes = PairingFacesBT[i][0]
    BSF = PairingFacesBT[i][1]
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-FaceBT-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-FaceTNode'+str(FaceTNodes[i]+1),j+1),
            (BSF[0],str(InstanceName)+'-FaceBNode'+str(BNodes[0]+1),j+1),
            (BSF[1],str(InstanceName)+'-FaceBNode'+str(BNodes[1]+1),j+1),
            (BSF[2],str(InstanceName)+'-FaceBNode'+str(BNodes[2]+1),j+1),
            (BSF[3],str(InstanceName)+'-FaceBNode'+str(BNodes[3]+1),j+1),
            (1.0,'V4',j+1),
            (-1.0,'V1',j+1)))

# FaceBa and FaceF
for i in range(len(FaceBaNodes)):
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceBaNode'+str(FaceBaNodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[FaceBaNodes[i]:FaceBaNodes[i]+1])
    
for i in range(len(FaceFNodes)):
    # Skip V1 and V5 which are control nodes
    if FaceFNodes[i] == V5:
        continue
    
    mdb.models[ModelName].rootAssembly.Set(name=str(InstanceName)+'-FaceFNode'+str(FaceFNodes[i]+1),nodes=mdb.models[ModelName].rootAssembly.instances[InstanceName].nodes[FaceFNodes[i]:FaceFNodes[i]+1])
    
    BaNodes = PairingFacesBaF[i][0]
    BaSF = PairingFacesBaF[i][1]
    
    for j in range(3):
        mdb.models[ModelName].Equation(name=str(InstanceName)+'-FaceBaF-'+str(i+1)+'-DOF-'+str(j+1),terms=(
            (-1.0,str(InstanceName)+'-FaceFNode'+str(FaceFNodes[i]+1),j+1),
            (BaSF[0],str(InstanceName)+'-FaceBaNode'+str(BaNodes[0]+1),j+1),
            (BaSF[1],str(InstanceName)+'-FaceBaNode'+str(BaNodes[1]+1),j+1),
            (BaSF[2],str(InstanceName)+'-FaceBaNode'+str(BaNodes[2]+1),j+1),
            (BaSF[3],str(InstanceName)+'-FaceBaNode'+str(BaNodes[3]+1),j+1),
            (1.0,'V5',j+1),
            (-1.0,'V1',j+1)))


'''
Proposed revisions (to be implemented)
Remove terms in the equation where coefficient is 0
Debug why loadings in all directions except uniaxial in x results in uneven deformation 
'''



