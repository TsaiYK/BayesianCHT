"""
Description:
This script solves the ABAQUS portion of Major Project. It takes in different
values for wing geometry and writes to a file in order to
view the output (max Mises stress, mass, eigenvalue for buckling, and identifying the design is successful or failed).

Script by:
Ying-Kuan Tsai, Lorenzo Casero, Haynes Herr, and Brendan White
Dept. of Aerospace Engineering
Texas A&M University
December 4th, 2021
"""

from abaqus import *
from abaqusConstants import *
import __main__
import section
import odbSection
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import displayGroupOdbToolset as dgo
import time
import numpy as np
from math import atan, sin, cos, tan
from EigenvalueExtraction import findEigenValue
# from DOEmethods import LHS
#from Post_P_Script import getResults

######################################
# Variable and Fixed Design Parameters
######################################

##########################
# FEA Modeling Parameters
# (e.g., mesh seeds, step times, etc)
##########################

session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)

####################################
### Calculated Properties/Values ###
####################################

filename = 'WingOutput.txt'
### Write data file column headings
DataFile = open(filename,'w')
# DataFile.write('numRibs, ribThick, SparThick, skinThick, stringThick, csf, MaterialSelection\n max Mises stress, mass, eig, buckle\n')
DataFile.close()

filenameDisppTipX = 'WingDisppTipX.txt'
### Write data file column headings
DataFile = open(filenameDisppTipX,'w')
DataFile.close()

filenameDisppTipY = 'WingDisppTipY.txt'
### Write data file column headings
DataFile = open(filenameDisppTipY,'w')
DataFile.close()

###############################
### Generation of FEA Model ###
###############################

### Note: If you create a loop, START it here 
### (i.e., the contents of the loop below are intended)
# seedSize = [0.25, 0.2, 0.15, 0.1, 0.075]
sizes = 0.02

### Write data file column headings

scaleFac = 1 # chord length (m)
# Load airfoil geometry by the x-y coordinates
x_coordinate = []
f = open('airfoil_NACA0012_20grids_x.txt')
for line in f.readlines():
	x_coordinate.append(float(line))
f.close()
numNodes = len(x_coordinate)

for i in range(numNodes-1):
	x_coordinate[i] = x_coordinate[i]*scaleFac
	
y_coordinate = []
f = open('airfoil_NACA0012_20grids_y.txt')
for line in f.readlines():
	y_coordinate.append(float(line))
f.close()
for i in range(numNodes-1):
	y_coordinate[i] = y_coordinate[i]*scaleFac

h = y_coordinate[5]-y_coordinate[13]
w = x_coordinate[5]-x_coordinate[13]

# calculate the slope of lines
slope = [0]*(numNodes-1)
for i in range(numNodes-1):
	slope[i] = atan((y_coordinate[i+1]-y_coordinate[i])/(x_coordinate[i+1]-x_coordinate[i]))*180/np.pi-90

entryPoint = 1
exitPoint = 19
horEntryPoint = 5
horExitPoint = 6
horEntryPoint_bot = 13
horExitPoint_bot = 14
points_airfoil1 = [[0]*2]*numNodes
k = 0
for i in range(entryPoint,horEntryPoint+1):
	points_airfoil1[k][0] = x_coordinate[i]
	points_airfoil1[k][1] = y_coordinate[i]
	points_airfoil1[k] = tuple(points_airfoil1[i])
	k = k+1
points_airfoil1[-1] = tuple(points_airfoil1[-1])

points_airfoil2 = [[0]*2]*numNodes
k = 0
for i in range(horExitPoint,horEntryPoint_bot+1):
	points_airfoil2[k][0] = x_coordinate[i]
	points_airfoil2[k][1] = y_coordinate[i]
	points_airfoil2[k] = tuple(points_airfoil2[i])
	k = k+1
points_airfoil2[-1] = tuple(points_airfoil2[-1])

points_airfoil3 = [[0]*2]*numNodes
k = 0
for i in range(horExitPoint_bot,exitPoint):
	points_airfoil3[k][0] = x_coordinate[i]
	points_airfoil3[k][1] = y_coordinate[i]
	points_airfoil3[k] = tuple(points_airfoil3[i])
	k = k+1
points_airfoil3[-1] = tuple(points_airfoil3[-1])

## Design Variables
# numRibs = 10 # (2 - 12)
# ribThick = 0.006143791649 # 0.001-0.01
# stringThick = 0.001039483248 # 0.0005-0.0015
# skinThick = 0.00127515534894539 # 0.0005-0.0015
# SparThick = 0.0195576853897116	# 0.01-0.02
# csf = 1.04424340902093	# clevis scale factor (1.0-1.4)
# MaterialSelection = 'Titanium'
# MaterialSelectionIndex = 1

selectExclude2 = 0.003
# Read the input data from the text file
nV = 7 # number of design variables
DV_file = 'DesignVariables.txt'

A = [0]*nV
i = 0
file_in = open(DV_file, 'r')
for y in file_in.read().split('\n'):
	if i<nV:
		A[i] = (float(y))
		i = i+1
file_in.close()
numRibs, ribThick, stringThick, skinThick, SparThick, csf, MaterialSelectionIndex = A
numRibs = int(numRibs)
MaterialSelectionIndex = int(MaterialSelectionIndex)

if MaterialSelectionIndex == 0:
	MaterialSelection = 'Aluminum'
	Density = 2780
	YoungsModulus = 73.1e9
	PoissonRatio = 0.33
	YieldTensile = 290000000
else:
	MaterialSelection = 'Titanium'
	Density = 4430
	YoungsModulus = 113.8e9
	PoissonRatio = 0.342
	YieldTensile = 828000000

### Independent
clevisSeed = 0.005
assemblySeed = 0.02
wingLength = 5.0


### Dependent DVs
# rib variables (radii)
t_rib = 0.01
t_rib_spar = 0.02
r = [0.015,0.015] # (m)
for i in range(len(r)):
	r[i] = r[i]*scaleFac
	
# skin variables
partLength = wingLength # wing length (m)

## Spar variables# Fixed:
SparWidth = w/2
CleviceLength = 0.15
SparHeight = h
SparLength = wingLength
selectExclude = CleviceLength*0.1

## Clevice variablesSparHeight = 60
Circlex = CleviceLength/5
Circley = SparHeight/7
Radius = 0.15*Circlex
filletRad = 0.003

## Stringer variables
length_string = wingLength
height_string = h/20
width_string = w/20

x_circle = [0.03,0.68]
for i in range(len(r)):
	x_circle[i] = x_circle[i]*scaleFac
	
y_circle = [0]*len(x_circle)
for i in range(len(r)):
	y_circle[i] = -0.0162*x_circle[i]

Mdb() 

tic = time.time()

nS = 1

for i in range(nS):
	print i

	### Scripting the entire model allows its entire
	### contents to be packaged into this single file. 
	  
	# Sketch Geometry and Create Parts
	print 'Sketching/Creating the part'
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=100.0)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=STANDALONE)
	
	## Airfoil rib
	# part 1
	for i in range(0,horEntryPoint+1):
		s.Line(point1=(x_coordinate[i], y_coordinate[i]), point2=(x_coordinate[i+1], y_coordinate[i+1]))
	# part 2
	s.Spline(points=tuple(points_airfoil2))
	# part 3: horizontal line
	s.Line(point1=(x_coordinate[horEntryPoint_bot], y_coordinate[horEntryPoint_bot]), point2=(x_coordinate[horExitPoint_bot], y_coordinate[horExitPoint_bot]))
	for i in range(horExitPoint_bot,exitPoint):
		s.Line(point1=(x_coordinate[i], y_coordinate[i]), point2=(x_coordinate[i+1], y_coordinate[i+1]))
	s.Line(point1=(x_coordinate[numNodes-1], y_coordinate[numNodes-1]), point2=(x_coordinate[0], y_coordinate[0]))
	
	# First hole (cicle)
	s.CircleByCenterPerimeter(center=(x_circle[0], y_circle[0]), point1=(x_circle[0], r[0]+y_circle[0]))
	# Second hole
	s.Line(point1=(0.07498, -0.04348+t_rib), point2=(0.07645, 0.04083-t_rib))
	s.Line(point1=(0.07645, 0.04083-t_rib), point2=(0.195-t_rib_spar, 0.05363-t_rib-t_rib_spar*0.1))
	s.Line(point1=(0.195-t_rib_spar, 0.05363-t_rib-t_rib_spar*0.1), point2=(0.193-t_rib_spar, -0.0604+t_rib+t_rib_spar*0.01))
	s.Line(point1=(0.193-t_rib_spar, -0.0604+t_rib+t_rib_spar*0.01), point2=(0.07498, -0.04348+t_rib))
	# Third hole
	s.Line(point1=(0.3349, -0.06559+t_rib), point2=(0.337, 0.05387-t_rib))
	s.Line(point1=(0.337, 0.05387-t_rib), point2=(0.6323, 0.03192-t_rib))
	s.Line(point1=(0.6323, 0.03192-t_rib), point2=(0.6308, -0.05397+t_rib))
	s.Line(point1=(0.6308, -0.05397+t_rib), point2=(0.3349, -0.06559+t_rib))
	# Forth hole (circle)
	s.CircleByCenterPerimeter(center=(x_circle[1], y_circle[1]), point1=(x_circle[1], y_circle[1]+r[1]))

	
	

	p = mdb.models['Model-1'].Part(name='Rib', dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	p = mdb.models['Model-1'].parts['Rib']
	p.BaseSolidExtrude(sketch=s, depth=0.01)
	s.unsetPrimaryObject()
	
	# Create shell from solid
	p = mdb.models['Model-1'].parts['Rib']
	c1 = p.cells
	p.RemoveCells(cellList=(c1.findAt(coordinates=(0.878811, -0.032036, 0.006667)), 
		))
	
	# Remove the front face
	p = mdb.models['Model-1'].parts['Rib']
	f = p.faces
	p.RemoveFaces(faceList=(f.findAt(coordinates=(0.662178, 0.012915, 0.01)), ), 
		deleteCells=False)
	
	
	# Remove the faces of holes
	p = mdb.models['Model-1'].parts['Rib']
	session.viewports['Viewport: 1'].setValues(displayedObject=p)
	session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
	
	p = mdb.models['Model-1'].parts['Rib']
	f1 = p.faces
	p.RemoveFaces(faceList=(f1.findAt(coordinates=(0.6313, -0.022007, 0.006667)), 
		f1.findAt(coordinates=(0.533867, 0.029237, 0.006667)), f1.findAt(
		coordinates=(0.3363, 0.010717, 0.006667)), f1.findAt(coordinates=(0.433533, 
		-0.051717, 0.006667)), f1.findAt(coordinates=(0.173667, -0.01959, 
		0.006667)), f1.findAt(coordinates=(0.14215, 0.03803, 0.006667)), f1.findAt(
		coordinates=(0.07596, 0.009393, 0.006667)), f1.findAt(coordinates=(
		0.107653, -0.039053, 0.006667)), f1.findAt(coordinates=(0.028053, 0.014387, 
		0.006667)), f1.findAt(coordinates=(0.678053, 0.003857, 0.006667))), 
		deleteCells=False)
	
	
	## Skin
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=100.0)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=STANDALONE)
	# Airfoil loop
	# part 1
	for i in range(0,horEntryPoint+1):
		s.Line(point1=(x_coordinate[i], y_coordinate[i]), point2=(x_coordinate[i+1], y_coordinate[i+1]))
	# part 2
	s.Spline(points=tuple(points_airfoil2))
	# part 3: horizontal line
	s.Line(point1=(x_coordinate[horEntryPoint_bot], y_coordinate[horEntryPoint_bot]), point2=(x_coordinate[horExitPoint_bot], y_coordinate[horExitPoint_bot]))
	for i in range(horExitPoint_bot,exitPoint):
		s.Line(point1=(x_coordinate[i], y_coordinate[i]), point2=(x_coordinate[i+1], y_coordinate[i+1]))
	s.Line(point1=(x_coordinate[numNodes-1], y_coordinate[numNodes-1]), point2=(x_coordinate[0], y_coordinate[0]))

	p = mdb.models['Model-1'].Part(name='Skin', dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	p = mdb.models['Model-1'].parts['Skin']
	# Airfoil Extruded
	p.BaseShellExtrude(sketch=s, depth=partLength)
	s.unsetPrimaryObject()
	# View and finish part
	session.viewports['Viewport: 1'].setValues(displayedObject=p)
	del mdb.models['Model-1'].sketches['__profile__']
	

	
	## Spar
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=100.0)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=STANDALONE)
	s.Line(point1=(0.0, 0.0), point2=(SparWidth, 0.0))
	s.HorizontalConstraint(entity=g.findAt((SparWidth/2, 0.0)), addUndoState=False)
	s.Line(point1=(SparWidth, 0.0), point2=(SparWidth, SparHeight))
	s.VerticalConstraint(entity=g.findAt((SparWidth, SparHeight/2)), addUndoState=False)
	s.PerpendicularConstraint(entity1=g.findAt((SparWidth/2, 0.0)), entity2=g.findAt((SparWidth, 
		SparHeight/2)), addUndoState=False)
	s.Line(point1=(SparWidth, SparHeight), point2=(0.0, SparHeight))
	s.HorizontalConstraint(entity=g.findAt((SparWidth/2, SparHeight)), addUndoState=False)
	s.PerpendicularConstraint(entity1=g.findAt((SparWidth, SparHeight/2)), entity2=g.findAt((
		SparWidth/2, SparHeight)), addUndoState=False)
	p = mdb.models['Model-1'].Part(name='Spar', dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	p = mdb.models['Model-1'].parts['Spar']
	p.BaseShellExtrude(sketch=s, depth=SparLength)
	s.unsetPrimaryObject()
	p = mdb.models['Model-1'].parts['Spar']
	
	
	# ## Clevice
	s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=200.0)
	g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
	s1.setPrimaryObject(option=STANDALONE)
	s1.Line(point1=(0.0, 0.0), point2=(0.0, SparHeight))
	s1.Line(point1=(0.0, 0.0), point2=(-SparWidth, 0.0))
	s1.Line(point1=(0.0, SparHeight), point2=(-SparWidth, SparHeight))
	s1.Line(point1=(-SparWidth, SparHeight), point2=(-SparWidth, SparHeight-SparThick))
	s1.Line(point1=(-SparWidth, SparHeight-SparThick), point2=(-SparThick, SparHeight-SparThick))
	s1.Line(point1=(-SparThick, SparHeight-SparThick), point2=(-SparThick, SparThick))
	s1.Line(point1=(-SparThick, SparThick), point2=(-SparWidth, SparThick))
	s1.Line(point1=(-SparWidth, SparThick), point2=(-SparWidth, 0.0))
	p = mdb.models['Model-1'].Part(name='Clevice', dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	p = mdb.models['Model-1'].parts['Clevice']
	p.BaseSolidExtrude(sketch=s1, depth=CleviceLength)
	s1.unsetPrimaryObject()
	p = mdb.models['Model-1'].parts['Clevice']
	
	# #Cut Extrudes:
	p = mdb.models['Model-1'].parts['Clevice']
	f, e1 = p.faces, p.edges
	t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(0.0, 0.5*h, 
	   0.5*CleviceLength)), sketchUpEdge=e1.findAt(coordinates=(0.0, 0.5*h, 0.0)), 
		sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, CleviceLength))
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=156.2, gridSpacing=3.9, transform=t)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=SUPERIMPOSE)
	p = mdb.models['Model-1'].parts['Clevice']
	p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

	#Holes:
	s.CircleByCenterPerimeter(center=(Circlex, Circley), point1=(Circlex+Radius, Circley))
	s.CircleByCenterPerimeter(center=(Circlex, SparHeight-Circley), point1=(Circlex+Radius,SparHeight-Circley))
	#Arcs:
	s.Arc3Points(point1=(Circlex, SparHeight), point2=(Circlex, SparHeight-2*Circley), point3=(Circlex-Circley, SparHeight-Circley))
	s.Arc3Points(point1=(Circlex, 0.0), point2=(Circlex, 2*Circley), point3=(Circlex-Circley, Circley))
	s.Arc3Points(point1=(Circlex, 2*Circley), point2=(Circlex, SparHeight-2*Circley), point3=(Circlex*2, SparHeight/2))	
	s.Line(point1=(Circlex, SparHeight), point2=(0.0, SparHeight))
	s.Line(point1=(0.0, SparHeight), point2=(0.0, 0.0))
	s.Line(point1=(0.0, 0.0), point2=(Circlex, 0.0))
	#Diamond:
	DX1 = CleviceLength*2/3*1/csf
	DY1 = SparHeight/2
	DX2 = (DX1+CleviceLength)/2
	DY2 = (3*SparHeight/4 - Circley)*csf
	DX3 = CleviceLength
	DY3 = (DY2+DY1)/2
	DX4 = DX3
	DY4 = SparHeight-(DY2+DY1)/2	
	DX5 = DX2
	DY5 = SparHeight-DY2
	DX6 = DX1
	DY6 = DY1
	s.Line(point1=(DX1, DY1), point2=(DX2, DY2))
	s.Line(point1=(DX2, DY2), point2=(DX3,DY3))
	s.Line(point1=(DX3, DY3), point2=(DX4, DY4))
	s.Line(point1=(DX4,DY4), point2=(DX5,DY5))
	s.Line(point1=(DX5, DY5), point2=(DX6, DY6))
	c1x = (DX1+DX2)/2
	c1y = (DY1+DY2)/2
	c2x = c1x
	c2y = SparHeight-(DY2+DY1)/2
	c3x = (CleviceLength+DX2)/2
	c3y = (DY2+(DY2+DY1)/2)/2
	c4x = (CleviceLength+DX2)/2
	c4y = (SparHeight-DY2+SparHeight-(DY1+DY2)/2)/2
	s.FilletByRadius(radius=filletRad, curve1=g.findAt((c1x,c1y)), 
		nearPoint1=(c1x, c1y), curve2=g.findAt((
		c2x, c2y)), nearPoint2=(c2x, c2y))
	s.FilletByRadius(radius=filletRad, curve1=g.findAt((c3x, c3y)), 
		nearPoint1=(c3x, c3y), curve2=g.findAt((
		c1x, c1y)), nearPoint2=(c1x, c1y))
	s.FilletByRadius(radius=filletRad, curve1=g.findAt((c2x, c2y)), 
		nearPoint1=(c2x, c2y), curve2=g.findAt((c4x, 
		c4y)), nearPoint2=(c4x, c4y))
	f, e = p.faces, p.edges
	p.CutExtrude(sketchPlane=f.findAt(coordinates=(0.0, 0.5*h, 
	   0.5*CleviceLength)), 
		sketchUpEdge=e.findAt(coordinates=(0.0, 0.5*h, 0.0)), sketchPlaneSide=SIDE1, 
		sketchOrientation=RIGHT, sketch=s, flipExtrudeDirection=OFF)
	# Partitions:
	p = mdb.models['Model-1'].parts['Clevice']
	c = p.cells
	pickedCells = c.findAt(((-0.071993, 0.107598, 0.085016), ))
	e1, v2, d2 = p.edges, p.vertices, p.datums
	p.PartitionCellByPlanePointNormal(point=v2.findAt(coordinates=(-SparThick, SparHeight-SparThick, 
		0.0)), normal=e1.findAt(coordinates=(0.0, Circley, 0.0)), 
		cells=pickedCells)
	p = mdb.models['Model-1'].parts['Clevice']
	c = p.cells
	pickedCells = c.findAt(((-0.006667, 0.090594, 0.0), ))
	e, v1, d1 = p.edges, p.vertices, p.datums
	p.PartitionCellByPlanePointNormal(point=v1.findAt(coordinates=(-SparThick, SparThick, 
		0.0)), normal=e.findAt(coordinates=(0.0, Circley, 0.0)), 
		cells=pickedCells)
	p = mdb.models['Model-1'].parts['Clevice']
	c = p.cells
	pickedCells = c.findAt(((0.0, 0.003333, 0.04), ), ((-0.003333, 0.10191, 
		0.11788), ), ((0.0, 0.109351, 0.128851), ))
	e1, v2, d2 = p.edges, p.vertices, p.datums
	p.PartitionCellByPlanePointNormal(normal=e1.findAt(coordinates=(0.0, SparHeight, 
		Circlex)), cells=pickedCells, point=p.InterestingPoint(edge=e1.findAt(
		coordinates=(0.0, 0.096403, 0.136251)), rule=CENTER))
	# Output Selection Partition/Set:
	p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=selectExclude)
	pickedCells = c[:]
	p.PartitionCellByDatumPlane(datumPlane=d1[6], cells=pickedCells)
	cells = c.getByBoundingBox(-1.0,-1.0,selectExclude-0.00001,1.0,1.0,1)
	region = p.Set(cells=cells, name='meshSet1')
	
	## Stringer
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=100.0)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=STANDALONE)
	s.Line(point1=(0.0, height_string), point2=(0.0, 0.0))
	s.Line(point1=(0.0, 0.0), point2=(width_string, 0.0))
	p = mdb.models['Model-1'].Part(name='Stringer', dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	p = mdb.models['Model-1'].parts['Stringer']
	p.BaseShellExtrude(sketch=s, depth=length_string)
	s.unsetPrimaryObject()
	p = mdb.models['Model-1'].parts['Stringer']
	del mdb.models['Model-1'].sketches['__profile__']

	#Defining the face partitions
	print 'Partitioning part'
	

	# Create Material
	print 'Creating the Materials'
	mdb.models['Model-1'].Material(name=MaterialSelection)
	mdb.models['Model-1'].materials[MaterialSelection].Density(table=((Density, ), ))
	mdb.models['Model-1'].materials[MaterialSelection].Elastic(table=((YoungsModulus, 
		PoissonRatio), ))
	
	#Create/Assign Section
	print 'Creating the Sections'
	print 'Assigning the Sections'
	
	# Rib
	p = mdb.models['Model-1'].parts['Rib']
	session.viewports['Viewport: 1'].setValues(displayedObject=p)
	mdb.models['Model-1'].HomogeneousShellSection(name='Rib-Sec', 
		preIntegrate=OFF, material=MaterialSelection, thicknessType=UNIFORM, 
		thickness=ribThick, thicknessField='', nodalThicknessField='', 
		idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
		thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
		integrationRule=SIMPSON, numIntPts=5)
	session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
	p = mdb.models['Model-1'].parts['Rib']
	f = p.faces
	faces = f.findAt(((0.97588, -0.021505, 0.006667), ), ((0.878811, -0.032036, 
		0.006667), ), ((0.730563, -0.04594, 0.006667), ), ((0.581118, -0.056781, 
		0.006667), ), ((0.432828, -0.061738, 0.006667), ), ((0.287576, -0.060397, 
		0.006667), ), ((0.158139, -0.056784, 0.006667), ), ((0.242298, 0.053867, 
		0.006667), ), ((0.385871, 0.051097, 0.006667), ), ((0.533217, 0.04101, 
		0.006667), ), ((0.681997, 0.026168, 0.006667), ), ((0.83042, 0.008005, 
		0.006667), ), ((0.952202, -0.00894, 0.006667), ), ((0.999855, -0.017032, 
		0.006667), ), ((0.666005, 0.012915, 0.0), ))
	region = p.Set(faces=faces, name='Set-2')
	p = mdb.models['Model-1'].parts['Rib']
	p.SectionAssignment(region=region, sectionName='Rib-Sec', offset=0.0, 
		offsetType=TOP_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)
		
	# Skin
	p = mdb.models['Model-1'].parts['Skin']
	session.viewports['Viewport: 1'].setValues(displayedObject=p)
	session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
	mdb.models['Model-1'].HomogeneousShellSection(name='Skin-Sec', 
		preIntegrate=OFF, material=MaterialSelection, thicknessType=UNIFORM, 
		thickness=skinThick, thicknessField='', nodalThicknessField='', 
		idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
		thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
		integrationRule=SIMPSON, numIntPts=5)
	p = mdb.models['Model-1'].parts['Skin']
	f = p.faces
	faces = f.findAt(((0.976036, -0.012566, 3.333333), ), ((0.879394, 0.001346, 
		3.333333), ), ((0.731721, 0.020416, 3.333333), ), ((0.582745, 0.036465, 
		3.333333), ), ((0.43478, 0.048326, 3.333333), ), ((0.289629, 0.053867, 
		3.333333), ), ((0.160031, 0.05138, 3.333333), ), ((0.240276, -0.060397, 
		3.333333), ), ((0.383852, -0.061067, 3.333333), ), ((0.531461, -0.059594, 
		3.333333), ), ((0.680669, -0.049954, 3.333333), ), ((0.829634, -0.036981, 
		3.333333), ), ((0.951934, -0.024297, 3.333333), ), ((0.99984, -0.017872, 
		3.333333), ))
	region = p.Set(faces=faces, name='Set-2')
	p = mdb.models['Model-1'].parts['Skin']
	p.SectionAssignment(region=region, sectionName='Skin-Sec', offset=0.0, 
		offsetType=BOTTOM_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)
	
	
			
	# Spar
	p = mdb.models['Model-1'].parts['Spar']
	mdb.models['Model-1'].HomogeneousShellSection(name='Spar-Sec', 
		preIntegrate=OFF, material=MaterialSelection, thicknessType=UNIFORM, 
		thickness=SparThick, thicknessField='', nodalThicknessField='', 
		idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
		thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
		integrationRule=SIMPSON, numIntPts=5)
	p = mdb.models['Model-1'].parts['Spar']
	f = p.faces
	faces = f.findAt(((0.023998, 0.114264, 3.333333), ), ((0.071993, 0.076176, 
		3.333333), ), ((0.047995, 0.0, 3.333333), ))
	region = p.Set(faces=faces, name='Set-2')
	p = mdb.models['Model-1'].parts['Spar']
	p.SectionAssignment(region=region, sectionName='Spar-Sec', offset=0.0, 
		offsetType=TOP_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)	 
	
	# Stringer
	p = mdb.models['Model-1'].parts['Stringer']
	mdb.models['Model-1'].HomogeneousShellSection(name='Stringer_sec', 
		preIntegrate=OFF, material=MaterialSelection, thicknessType=UNIFORM, 
		thickness=stringThick, thicknessField='', nodalThicknessField='', 
		idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
		thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
		integrationRule=SIMPSON, numIntPts=5)
	p = mdb.models['Model-1'].parts['Stringer']
	f = p.faces
	faces = f.findAt(((0.0048, 0.0, 3.333333), ), ((0.0, 0.001904, 3.333333), ))
	region = p.Set(faces=faces, name='Set-2')
	p = mdb.models['Model-1'].parts['Stringer']
	p.SectionAssignment(region=region, sectionName='Rib-Sec', offset=0.0, 
		offsetType=TOP_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)
		
	# Clevice
	p = mdb.models['Model-1'].parts['Clevice']
	mdb.models['Model-1'].HomogeneousSolidSection(name='Clevice-Sec', 
		material=MaterialSelection, thickness=None)
	p = mdb.models['Model-1'].parts['Clevice']
	c = p.cells
	# cells = c.findAt(((-0.071993, 0.110521, 0.122074), ), ((-0.006667, 0.093463, 
		# 0.12044), ), ((0.0, 0.110931, 0.08), ), ((-0.006667, 0.035233, 0.04887), 
		# ), ((-0.006667, 0.011845, 0.12044), ), ((0.0, 0.003744, 0.122074), ), ((
		# -0.023998, 0.0, 0.08), ))
	cells = c.getByBoundingBox(-1.0,-1.0,-1.0,1.0,1.0,1.0)
	region = p.Set(cells=cells, name='Set-3')
	p = mdb.models['Model-1'].parts['Clevice']
	p.SectionAssignment(region=region, sectionName='Clevice-Sec', offset=0.0, 
		offsetType=MIDDLE_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)	

	# Assemble Parts
	print 'Placing Parts in Space'
	#Create Instances here
	# Rib and skin
	a1 = mdb.models['Model-1'].rootAssembly
	p = mdb.models['Model-1'].parts['Skin']
	a1.Instance(name='Skin-1', part=p, dependent=ON)
	
	for i in range(numRibs):
		RibName = 'Rib-'+str(i)
		a1 = mdb.models['Model-1'].rootAssembly
		p = mdb.models['Model-1'].parts['Rib']
		a1.Instance(name=RibName, part=p, dependent=ON)
		a1.translate(instanceList=(RibName, ), vector=(0.0, 0.0, wingLength/(numRibs-1)*i))
	
	
	# Spar
	a1 = mdb.models['Model-1'].rootAssembly
	a = mdb.models['Model-1'].rootAssembly
	a1 = mdb.models['Model-1'].rootAssembly
	p = mdb.models['Model-1'].parts['Spar']
	a1.Instance(name='Spar-1', part=p, dependent=ON)
	a1 = mdb.models['Model-1'].rootAssembly
	a1.rotate(instanceList=('Spar-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(
		0.0, 1.0, 0.0), angle=180.0)
	#: The instance Spar-1 was rotated by 180. degrees about the axis defined by the point 0., 0., 0. and the vector 0., 1., 0.
	a1 = mdb.models['Model-1'].rootAssembly
	a1.translate(instanceList=('Spar-1', ), vector=(0.1929+SparWidth, 0.0, 0.0))
	#: The instance Spar-1 was translated by 500.E-03, 0., 0. with respect to the assembly coordinate system
	a1 = mdb.models['Model-1'].rootAssembly
	a1.translate(instanceList=('Spar-1', ), vector=(0.0, -0.0604, 0.0))
	#: The instance Spar-1 was translated by 0., -200.E-03, 0. with respect to the assembly coordinate system
	
	# Clevice
	a = mdb.models['Model-1'].rootAssembly
	a = mdb.models['Model-1'].rootAssembly
	p = mdb.models['Model-1'].parts['Clevice']
	a.Instance(name='Clevice-1', part=p, dependent=ON)
	a = mdb.models['Model-1'].rootAssembly
	a.rotate(instanceList=('Clevice-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(
		0.0, 1.0, 0.0), angle=180.0)
	#: The instance Spar-1 was rotated by 180. degrees about the axis defined by the point 0., 0., 0. and the vector 0., 1., 0.
	a = mdb.models['Model-1'].rootAssembly
	a.translate(instanceList=('Clevice-1', ), vector=(0.1929, 0.0, 0.0))
	#: The instance Spar-1 was translated by 500.E-03, 0., 0. with respect to the assembly coordinate system
	a = mdb.models['Model-1'].rootAssembly
	a.translate(instanceList=('Clevice-1', ), vector=(0.0, -0.0604, 0.0))
	#: The instance Spar-1 was translated by 0., -200.E-03, 0. with respect to the assembly coordinate system

	p = mdb.models['Model-1'].parts['Stringer']
	a1 = mdb.models['Model-1'].rootAssembly
	a1.regenerate()
	a = mdb.models['Model-1'].rootAssembly
	a1 = mdb.models['Model-1'].rootAssembly
	a1.translate(instanceList=('Spar-1', ), vector=(0.0, 0.0, wingLength))

	# Stringer
	a1 = mdb.models['Model-1'].rootAssembly
	p = mdb.models['Model-1'].parts['Stringer']
	a1.Instance(name='Stringer-1', part=p, dependent=ON)
	a1 = mdb.models['Model-1'].rootAssembly
	# a1.rotate(instanceList=('Stringer-1', ), axisPoint=(0.0, 0.0, 0.0), 
		# axisDirection=(0.0, 0.0, 1.0), angle=slope[0])
	# a1.translate(instanceList=('Stringer-1', ), vector=(x_coordinate[1], y_coordinate[1], 0.0))
	# a = mdb.models['Model-1'].rootAssembly
	for i in range(5):
		StringerName = 'Stringer-'+str(i)
		a1 = mdb.models['Model-1'].rootAssembly
		p = mdb.models['Model-1'].parts['Stringer']
		a1.Instance(name=StringerName, part=p, dependent=ON)
		a1.rotate(instanceList=(StringerName, ), axisPoint=(0.0, 0.0, 0.0), 
			axisDirection=(0.0, 0.0, 1.0), angle=slope[i])
		a1.translate(instanceList=(StringerName, ), vector=(x_coordinate[i+1], y_coordinate[i+1], 0.0))
	
	for i in range(5):
		StringerName = 'Stringer-'+str(i+5)
		a1 = mdb.models['Model-1'].rootAssembly
		p = mdb.models['Model-1'].parts['Stringer']
		a1.Instance(name=StringerName, part=p, dependent=ON)
		a1.rotate(instanceList=(StringerName, ), axisPoint=(0.0, 0.0, 0.0), 
			axisDirection=(0.0, 0.0, 1.0), angle=90+slope[18-i])
		a1.translate(instanceList=(StringerName, ), vector=(x_coordinate[18-i], y_coordinate[18-i], 0.0))
	
	# Define Sets
	print 'Defining Sets'	
	# for-loop picking points for set
	a = mdb.models['Model-1'].rootAssembly
	zPosition = [0]*numRibs
	for i in range(numRibs):
		RibName = 'Rib-'+str(i)
		RibNameSet = 'Rib-OML-'+str(i)
		v1 = a.instances[RibName].vertices
		zPosition[i] = wingLength/(numRibs-1)*i
		verts1 = v1.findAt(((0.999826, -0.018712, zPosition[i]), ), ((0.927989, -0.02709, 
			zPosition[i]), ), ((0.780457, -0.041927, zPosition[i]), ), ((0.630774, -0.053967, zPosition[i]), ), 
			((0.481804, -0.062408, zPosition[i]), ), ((0.334876, -0.060397, zPosition[i]), ), ((
			0.192975, -0.060397, zPosition[i]), ), ((0.194966, 0.053867, zPosition[i]), ), ((0.336961, 
			0.053867, zPosition[i]), ), ((0.483689, 0.045555, zPosition[i]), ), ((0.632273, 0.03192, 
			zPosition[i]), ), ((0.781445, 0.014664, zPosition[i]), ), ((0.928369, -0.005313, zPosition[i]), ), 
			((0.99987, -0.016193, zPosition[i]), ))
		a.Set(vertices=verts1, name=RibNameSet)
	
	# Tip node
	a = mdb.models['Model-1'].rootAssembly
	v1 = a.instances['Spar-1'].vertices
	verts1 = v1.findAt(((0.264893, 0.053864, 5.0), ))
	a.Set(vertices=verts1, name='TIPNODE')
	
	#Define BCs
	print 'Defining all BCs'
	a = mdb.models['Model-1'].rootAssembly
	f1 = a.instances['Clevice-1'].faces
	faces1 = f1.findAt(((0.196233, 0.04151, -0.11788), ), ((0.196233, -0.040107, 
		-0.11788), ))
	region = a.Set(faces=faces1, name='Set-1')
	mdb.models['Model-1'].EncastreBC(name='ClevisBC', createStepName='Initial', 
		region=region, localCsys=None)
	
	
	print 'Merge assemble wing'
	a = mdb.models['Model-1'].rootAssembly
	if numRibs == 2:

		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(
			a.instances['Skin-1'], a.instances['Rib-1'], 
			a.instances['Rib-0'], a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, 
			originalInstances=SUPPRESS, domain=GEOMETRY)

	elif numRibs == 3:

		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(
			a.instances['Skin-1'], a.instances['Rib-2'], a.instances['Rib-1'], 
			a.instances['Rib-0'], a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, 
			originalInstances=SUPPRESS, domain=GEOMETRY)

	elif numRibs == 4:

		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(
			a.instances['Skin-1'],	a.instances['Rib-2'], a.instances['Rib-1'], 
			a.instances['Rib-0'], a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, 
			originalInstances=SUPPRESS, domain=GEOMETRY)

	elif numRibs == 5:

		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(
			a.instances['Skin-1'],	a.instances['Rib-4'], 
			a.instances['Rib-3'], a.instances['Rib-2'], a.instances['Rib-1'], 
			a.instances['Rib-0'], a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, 
			originalInstances=SUPPRESS, domain=GEOMETRY)

	elif numRibs == 6:

		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(
			a.instances['Skin-1'], a.instances['Rib-5'], a.instances['Rib-4'], 
			a.instances['Rib-3'], a.instances['Rib-2'], a.instances['Rib-1'], 
			a.instances['Rib-0'], a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, 
			originalInstances=SUPPRESS, domain=GEOMETRY)

	elif numRibs == 7:
		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(a.instances['Skin-1'], 
			a.instances['Rib-0'], a.instances['Rib-1'], a.instances['Rib-2'], 
			a.instances['Rib-3'], a.instances['Rib-4'], a.instances['Rib-5'], 
			a.instances['Rib-6'], 
			a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
	elif numRibs == 8:
		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(a.instances['Skin-1'], 
			a.instances['Rib-0'], a.instances['Rib-1'], a.instances['Rib-2'], 
			a.instances['Rib-3'], a.instances['Rib-4'], a.instances['Rib-5'], 
			a.instances['Rib-6'], a.instances['Rib-7'], 
			a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
	elif numRibs == 9:
		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(a.instances['Skin-1'], 
			a.instances['Rib-0'], a.instances['Rib-1'], a.instances['Rib-2'], 
			a.instances['Rib-3'], a.instances['Rib-4'], a.instances['Rib-5'], 
			a.instances['Rib-6'], a.instances['Rib-7'], a.instances['Rib-8'], 
			a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
	elif numRibs == 10:
		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(a.instances['Skin-1'], 
			a.instances['Rib-0'], a.instances['Rib-1'], a.instances['Rib-2'], 
			a.instances['Rib-3'], a.instances['Rib-4'], a.instances['Rib-5'], 
			a.instances['Rib-6'], a.instances['Rib-7'], a.instances['Rib-8'], a.instances['Rib-9'], 
			a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
	elif numRibs == 11:
		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(a.instances['Skin-1'], 
			a.instances['Rib-0'], a.instances['Rib-1'], a.instances['Rib-2'], 
			a.instances['Rib-3'], a.instances['Rib-4'], a.instances['Rib-5'], 
			a.instances['Rib-6'], a.instances['Rib-7'], a.instances['Rib-8'], 
			a.instances['Rib-9'], a.instances['Rib-10'], 
			a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
	elif numRibs == 12:
		a.InstanceFromBooleanMerge(name='assemble_wing', instances=(a.instances['Skin-1'], 
			a.instances['Rib-0'], a.instances['Rib-1'], a.instances['Rib-2'], 
			a.instances['Rib-3'], a.instances['Rib-4'], a.instances['Rib-5'], 
			a.instances['Rib-6'], a.instances['Rib-7'], a.instances['Rib-8'], 
			a.instances['Rib-9'], a.instances['Rib-9'], a.instances['Rib-9'], 
			a.instances['Rib-10'], a.instances['Rib-11'], 
			a.instances['Spar-1'], a.instances['Stringer-0'], 
			a.instances['Stringer-1'], a.instances['Stringer-2'], 
			a.instances['Stringer-3'], a.instances['Stringer-4'], 
			a.instances['Stringer-5'], a.instances['Stringer-6'], 
			a.instances['Stringer-7'], a.instances['Stringer-8'], 
			a.instances['Stringer-9'], ), keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)


	
	#Define Steps
	print 'Defining Steps'
	mdb.models['Model-1'].BuckleStep(name='Buckle-Step', previous='Initial', 
		numEigen=1, eigensolver=LANCZOS, minEigen=0.0, blockSize=DEFAULT, 
		maxBlocks=DEFAULT)
	mdb.models['Model-1'].StaticStep(name='Load-Step', previous='Initial')

		
	#Create Loads
	print 'Defining Loads'
	# Generate Pressures
	Pressures = []
	f = open('Pressures.txt')
	for line in f.readlines():
		Pressures.append(float(line))
	f.close()

	TPressures = []
	for i in range(0,len(Pressures)-1):
		TPressures.append((Pressures[i]+Pressures[i+1])/2)
	TPressures.append((Pressures[19]+Pressures[0])/2)	 

	TTP=[]
	for i in range(0,6):
		TTP.append(TPressures[i])
	TTP.append((TPressures[6]+TPressures[7]+TPressures[8]+TPressures[9]+TPressures[10]+TPressures[11]+TPressures[12])/7)
	for i in range(13,20):
		TTP.append(TPressures[i])

	# Pressure 1
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.930251, -0.0056, 5.003333), ))
	side2Faces1 = s1.findAt(((0.932134, -0.005886, 0.006667), ), ((0.932134, 
		-0.005886, 0.67), ), ((0.932134, -0.005886, 1.006667), ), ((0.932134, 
		-0.005886, 1.67), ), ((0.932134, -0.005886, 2.006667), ), ((0.932134, 
		-0.005886, 2.67), ), ((0.932134, -0.005886, 3.006667), ), ((0.932134, 
		-0.005886, 3.67), ), ((0.932134, -0.005886, 4.006667), ), ((0.932134, 
		-0.005886, 4.67), ), ((0.955968, -0.009513, 0.34), ), ((0.955968, 
		-0.009513, 1.34), ), ((0.955968, -0.009513, 1.003333), ), ((0.955968, 
		-0.009513, 2.34), ), ((0.955968, -0.009513, 2.003333), ), ((0.955968, 
		-0.009513, 3.34), ), ((0.955968, -0.009513, 3.003333), ), ((0.955968, 
		-0.009513, 4.34), ), ((0.955968, -0.009513, 4.003333), ), ((0.955968, 
		-0.009513, 0.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-1')
	mdb.models['Model-1'].Pressure(name='Load-1', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[0], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-1', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[0], 
		amplitude=UNSET)		
	# Pressure 2
	side1Faces1 = s1.findAt(((0.783332, 0.014408, 5.003333), ))
	side2Faces1 = s1.findAt(((0.785219, 0.014151, 0.006667), ), ((0.785219, 
		0.014151, 0.67), ), ((0.785219, 0.014151, 1.006667), ), ((0.785219, 
		0.014151, 1.67), ), ((0.785219, 0.014151, 2.006667), ), ((0.785219, 
		0.014151, 2.67), ), ((0.785219, 0.014151, 3.006667), ), ((0.785219, 
		0.014151, 3.67), ), ((0.785219, 0.014151, 4.006667), ), ((0.785219, 
		0.014151, 4.67), ), ((0.834194, 0.007492, 0.34), ), ((0.834194, 0.007492, 
		1.34), ), ((0.834194, 0.007492, 1.003333), ), ((0.834194, 0.007492, 2.34), 
		), ((0.834194, 0.007492, 2.003333), ), ((0.834194, 0.007492, 3.34), ), ((
		0.834194, 0.007492, 3.003333), ), ((0.834194, 0.007492, 4.34), ), ((
		0.834194, 0.007492, 4.003333), ), ((0.834194, 0.007492, 0.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-2')
	mdb.models['Model-1'].Pressure(name='Load-2', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[1], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-2', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[1], 
		amplitude=UNSET)
	# Pressure 3
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.634165, 0.031701, 5.003333), ))
	side2Faces1 = s1.findAt(((0.636057, 0.031482, 0.006667), ), ((0.636057, 
		0.031482, 0.67), ), ((0.636057, 0.031482, 1.006667), ), ((0.636057, 
		0.031482, 1.67), ), ((0.636057, 0.031482, 2.006667), ), ((0.636057, 
		0.031482, 2.67), ), ((0.636057, 0.031482, 3.006667), ), ((0.636057, 
		0.031482, 3.67), ), ((0.636057, 0.031482, 4.006667), ), ((0.636057, 
		0.031482, 4.67), ), ((0.685781, 0.02573, 0.34), ), ((0.685781, 0.02573, 
		1.34), ), ((0.685781, 0.02573, 1.003333), ), ((0.685781, 0.02573, 2.34), ), 
		((0.685781, 0.02573, 2.003333), ), ((0.685781, 0.02573, 3.34), ), ((
		0.685781, 0.02573, 3.003333), ), ((0.685781, 0.02573, 4.34), ), ((0.685781, 
		0.02573, 4.003333), ), ((0.685781, 0.02573, 0.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-3')
	mdb.models['Model-1'].Pressure(name='Load-3', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[2], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-3', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[2], 
		amplitude=UNSET)
	# Pressure 4
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.485585, 0.045381, 5.003333), ))
	side2Faces1 = s1.findAt(((0.487482, 0.045207, 0.006667), ), ((0.487482, 
		0.045207, 0.67), ), ((0.487482, 0.045207, 1.006667), ), ((0.487482, 
		0.045207, 1.67), ), ((0.487482, 0.045207, 2.006667), ), ((0.487482, 
		0.045207, 2.67), ), ((0.487482, 0.045207, 3.006667), ), ((0.487482, 
		0.045207, 3.67), ), ((0.487482, 0.045207, 4.006667), ), ((0.487482, 
		0.045207, 4.67), ), ((0.53701, 0.040662, 0.34), ), ((0.53701, 0.040662, 
		1.34), ), ((0.53701, 0.040662, 1.003333), ), ((0.53701, 0.040662, 2.34), ), 
		((0.53701, 0.040662, 2.003333), ), ((0.53701, 0.040662, 3.34), ), ((
		0.53701, 0.040662, 3.003333), ), ((0.53701, 0.040662, 4.34), ), ((0.53701, 
		0.040662, 4.003333), ), ((0.53701, 0.040662, 0.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-4')
	mdb.models['Model-1'].Pressure(name='Load-4', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[3], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-4', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[3], 
		amplitude=UNSET)
	# Pressure 5
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.338863, 0.05376, 5.003333), ))
	side2Faces1 = s1.findAt(((0.340764, 0.053652, 0.006667), ), ((0.340764, 
		0.053652, 0.67), ), ((0.340764, 0.053652, 1.006667), ), ((0.340764, 
		0.053652, 1.67), ), ((0.340764, 0.053652, 2.006667), ), ((0.340764, 
		0.053652, 2.67), ), ((0.340764, 0.053652, 3.006667), ), ((0.340764, 
		0.053652, 3.67), ), ((0.340764, 0.053652, 4.006667), ), ((0.340764, 
		0.053652, 4.67), ), ((0.389673, 0.050881, 0.34), ), ((0.389673, 0.050881, 
		1.34), ), ((0.389673, 0.050881, 1.003333), ), ((0.389673, 0.050881, 2.34), 
		), ((0.389673, 0.050881, 2.003333), ), ((0.389673, 0.050881, 3.34), ), ((
		0.389673, 0.050881, 3.003333), ), ((0.389673, 0.050881, 4.34), ), ((
		0.389673, 0.050881, 4.003333), ), ((0.389673, 0.050881, 0.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-7')
	mdb.models['Model-1'].Pressure(name='Load-5', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[4], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-5', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[4], 
		amplitude=UNSET)
	# Pressure 6
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.242298, 0.053867, 5.006667), ))
	side2Faces1 = s1.findAt(((0.242298, 0.053867, 0.34), ), ((0.242298, 0.053867, 
		1.34), ), ((0.242298, 0.053867, 1.003333), ), ((0.242298, 0.053867, 2.34), 
		), ((0.242298, 0.053867, 2.003333), ), ((0.242298, 0.053867, 3.34), ), ((
		0.242298, 0.053867, 3.003333), ), ((0.242298, 0.053867, 4.34), ), ((
		0.242298, 0.053867, 4.003333), ), ((0.289629, 0.053867, 0.006667), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-8')
	mdb.models['Model-1'].Pressure(name='Load-6', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[5], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-6', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[5], 
		amplitude=UNSET)
	# Pressure 7
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.022632, -0.025555, 5.006667), ))
	side2Faces1 = s1.findAt(((0.192925, -0.060392, 0.006667), ), ((0.158159, 
		0.051232, 0.006667), ), ((0.194248, 0.053819, 0.006667), ), ((0.192925, 
		-0.060392, 0.67), ), ((0.158159, 0.051232, 0.67), ), ((0.194248, 0.053819, 
		0.67), ), ((0.192925, -0.060392, 1.006667), ), ((0.158159, 0.051232, 
		1.006667), ), ((0.194248, 0.053819, 1.006667), ), ((0.192925, -0.060392, 
		1.67), ), ((0.158159, 0.051232, 1.67), ), ((0.194248, 0.053819, 1.67), ), (
		(0.192925, -0.060392, 2.006667), ), ((0.158159, 0.051232, 2.006667), ), ((
		0.194248, 0.053819, 2.006667), ), ((0.192925, -0.060392, 2.67), ), ((
		0.158159, 0.051232, 2.67), ), ((0.194248, 0.053819, 2.67), ), ((0.192925, 
		-0.060392, 3.006667), ), ((0.158159, 0.051232, 3.006667), ), ((0.194248, 
		0.053819, 3.006667), ), ((0.192925, -0.060392, 3.67), ), ((0.158159, 
		0.051232, 3.67), ), ((0.194248, 0.053819, 3.67), ), ((0.192925, -0.060392, 
		4.006667), ), ((0.158159, 0.051232, 4.006667), ), ((0.194248, 0.053819, 
		4.006667), ), ((0.192925, -0.060392, 4.67), ), ((0.158159, 0.051232, 4.67), 
		), ((0.194248, 0.053819, 4.67), ), ((0.194937, 0.053865, 0.003333), ), ((
		0.194937, 0.053865, 1.34), ), ((0.194937, 0.053865, 0.34), ), ((0.194937, 
		0.053865, 2.34), ), ((0.194937, 0.053865, 1.003333), ), ((0.194937, 
		0.053865, 3.34), ), ((0.194937, 0.053865, 2.003333), ), ((0.194937, 
		0.053865, 4.34), ), ((0.194937, 0.053865, 3.003333), ), ((0.194937, 
		0.053865, 4.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-9')
	mdb.models['Model-1'].Pressure(name='Load-7', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[6], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-7', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[6], 
		amplitude=UNSET)
	# Pressure 8
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.287576, -0.060397, 5.006667), ))
	side2Faces1 = s1.findAt(((0.240276, -0.060397, 0.006667), ), ((0.287576, 
		-0.060397, 1.34), ), ((0.287576, -0.060397, 0.34), ), ((0.287576, 
		-0.060397, 2.34), ), ((0.287576, -0.060397, 1.003333), ), ((0.287576, 
		-0.060397, 3.34), ), ((0.287576, -0.060397, 2.003333), ), ((0.287576, 
		-0.060397, 4.34), ), ((0.287576, -0.060397, 3.003333), ), ((0.287576, 
		-0.060397, 4.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-10')
	mdb.models['Model-1'].Pressure(name='Load-8', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[7], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-8', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[7], 
		amplitude=UNSET)
	# Pressure 9
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.435228, -0.061771, 5.003333), ))
	side2Faces1 = s1.findAt(((0.388651, -0.061133, 0.006667), ), ((0.388651, 
		-0.061133, 0.67), ), ((0.388651, -0.061133, 1.006667), ), ((0.388651, 
		-0.061133, 1.67), ), ((0.388651, -0.061133, 2.006667), ), ((0.388651, 
		-0.061133, 2.67), ), ((0.388651, -0.061133, 3.006667), ), ((0.388651, 
		-0.061133, 3.67), ), ((0.388651, -0.061133, 4.006667), ), ((0.388651, 
		-0.061133, 4.67), ), ((0.339675, -0.060463, 0.003333), ), ((0.339675, 
		-0.060463, 1.34), ), ((0.339675, -0.060463, 0.34), ), ((0.339675, 
		-0.060463, 2.34), ), ((0.339675, -0.060463, 1.003333), ), ((0.339675, 
		-0.060463, 3.34), ), ((0.339675, -0.060463, 2.003333), ), ((0.339675, 
		-0.060463, 4.34), ), ((0.339675, -0.060463, 3.003333), ), ((0.339675, 
		-0.060463, 4.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-11')
	mdb.models['Model-1'].Pressure(name='Load-9', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[8], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-9', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[8], 
		amplitude=UNSET)
	# Pressure 10
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.583513, -0.056645, 5.003333), ))
	side2Faces1 = s1.findAt(((0.536253, -0.059323, 0.006667), ), ((0.536253, 
		-0.059323, 0.67), ), ((0.536253, -0.059323, 1.006667), ), ((0.536253, 
		-0.059323, 1.67), ), ((0.536253, -0.059323, 2.006667), ), ((0.536253, 
		-0.059323, 2.67), ), ((0.536253, -0.059323, 3.006667), ), ((0.536253, 
		-0.059323, 3.67), ), ((0.536253, -0.059323, 4.006667), ), ((0.536253, 
		-0.059323, 4.67), ), ((0.486596, -0.062137, 0.003333), ), ((0.486596, 
		-0.062137, 1.34), ), ((0.486596, -0.062137, 0.34), ), ((0.486596, 
		-0.062137, 2.34), ), ((0.486596, -0.062137, 1.003333), ), ((0.486596, 
		-0.062137, 3.34), ), ((0.486596, -0.062137, 2.003333), ), ((0.486596, 
		-0.062137, 4.34), ), ((0.486596, -0.062137, 3.003333), ), ((0.486596, 
		-0.062137, 4.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-12')
	mdb.models['Model-1'].Pressure(name='Load-10', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[9], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-10', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[9], 
		amplitude=UNSET)
	# Pressure 11
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.732955, -0.045748, 5.003333), ))
	side2Faces1 = s1.findAt(((0.685453, -0.049569, 0.006667), ), ((0.685453, 
		-0.049569, 0.67), ), ((0.685453, -0.049569, 1.006667), ), ((0.685453, 
		-0.049569, 1.67), ), ((0.685453, -0.049569, 2.006667), ), ((0.685453, 
		-0.049569, 2.67), ), ((0.685453, -0.049569, 3.006667), ), ((0.685453, 
		-0.049569, 3.67), ), ((0.685453, -0.049569, 4.006667), ), ((0.685453, 
		-0.049569, 4.67), ), ((0.635558, -0.053582, 0.003333), ), ((0.635558, 
		-0.053582, 1.34), ), ((0.635558, -0.053582, 0.34), ), ((0.635558, 
		-0.053582, 2.34), ), ((0.635558, -0.053582, 1.003333), ), ((0.635558, 
		-0.053582, 3.34), ), ((0.635558, -0.053582, 2.003333), ), ((0.635558, 
		-0.053582, 4.34), ), ((0.635558, -0.053582, 3.003333), ), ((0.635558, 
		-0.053582, 4.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-13')
	mdb.models['Model-1'].Pressure(name='Load-11', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[10], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-11', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[10], 
		amplitude=UNSET)
	# Pressure 12
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.881199, -0.031795, 5.003333), ))
	side2Faces1 = s1.findAt(((0.83441, -0.036501, 0.006667), ), ((0.83441, 
		-0.036501, 0.67), ), ((0.83441, -0.036501, 1.006667), ), ((0.83441, 
		-0.036501, 1.67), ), ((0.83441, -0.036501, 2.006667), ), ((0.83441, 
		-0.036501, 2.67), ), ((0.83441, -0.036501, 3.006667), ), ((0.83441, 
		-0.036501, 3.67), ), ((0.83441, -0.036501, 4.006667), ), ((0.83441, 
		-0.036501, 4.67), ), ((0.785233, -0.041447, 0.003333), ), ((0.785233, 
		-0.041447, 1.34), ), ((0.785233, -0.041447, 0.34), ), ((0.785233, 
		-0.041447, 2.34), ), ((0.785233, -0.041447, 1.003333), ), ((0.785233, 
		-0.041447, 3.34), ), ((0.785233, -0.041447, 2.003333), ), ((0.785233, 
		-0.041447, 4.34), ), ((0.785233, -0.041447, 3.003333), ), ((0.785233, 
		-0.041447, 4.003333), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-14')
	mdb.models['Model-1'].Pressure(name='Load-12', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[11], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-12', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[11], 
		amplitude=UNSET)
	# Pressure 13
	s1 = a.instances['assemble_wing-1'].faces
	side1Faces1 = s1.findAt(((0.978264, -0.021227, 5.003333), ))
	side2Faces1 = s1.findAt(((0.932756, -0.026534, 4.34), ), ((0.932756, -0.026534, 
		4.003333), ), ((0.932756, -0.026534, 3.34), ), ((0.932756, -0.026534, 
		3.003333), ), ((0.932756, -0.026534, 2.34), ), ((0.932756, -0.026534, 
		2.003333), ), ((0.932756, -0.026534, 1.34), ), ((0.932756, -0.026534, 
		1.003333), ), ((0.932756, -0.026534, 0.34), ), ((0.932756, -0.026534, 
		0.003333), ), ((0.956702, -0.023741, 0.67), ), ((0.956702, -0.023741, 
		1.67), ), ((0.956702, -0.023741, 1.006667), ), ((0.956702, -0.023741, 
		2.67), ), ((0.956702, -0.023741, 2.006667), ), ((0.956702, -0.023741, 
		3.67), ), ((0.956702, -0.023741, 3.006667), ), ((0.956702, -0.023741, 
		4.67), ), ((0.956702, -0.023741, 4.006667), ), ((0.956702, -0.023741, 
		0.006667), ))
	region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
		name='Surf-15')
	mdb.models['Model-1'].Pressure(name='Load-13', createStepName='Load-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[12], 
		amplitude=UNSET)
	mdb.models['Model-1'].Pressure(name='BLoad-13', createStepName='Buckle-Step', 
		region=region, distributionType=UNIFORM, field='', magnitude=TTP[12], 
		amplitude=UNSET)
	# Pressure 14
	if TTP[13] != 0:
		s1 = a.instances['assemble_wing-1'].faces
		side1Faces1 = s1.findAt(((0.999855, -0.017032, 5.006667), ))
		side2Faces1 = s1.findAt(((0.999855, -0.017032, 0.34), ), ((0.999855, -0.017032, 
			1.34), ), ((0.999855, -0.017032, 1.003333), ), ((0.999855, -0.017032, 
			2.34), ), ((0.999855, -0.017032, 2.003333), ), ((0.999855, -0.017032, 
			3.34), ), ((0.999855, -0.017032, 3.003333), ), ((0.999855, -0.017032, 
			4.34), ), ((0.999855, -0.017032, 4.003333), ), ((0.99984, -0.017872, 
			0.006667), ))
		region = a.Surface(side1Faces=side1Faces1, side2Faces=side2Faces1, 
			name='Surf-16')
		mdb.models['Model-1'].Pressure(name='Load-14', createStepName='Load-Step', 
			region=region, distributionType=UNIFORM, field='', magnitude=TTP[13], 
			amplitude=UNSET)
		mdb.models['Model-1'].Pressure(name='BLoad-14', createStepName='Buckle-Step', 
			region=region, distributionType=UNIFORM, field='', magnitude=TTP[13], 
			amplitude=UNSET)
	
	print 'Define Interaction and Add Constraints'
	
	a = mdb.models['Model-1'].rootAssembly
	a.regenerate()
	p = mdb.models['Model-1'].parts['Clevice']
	p = mdb.models['Model-1'].parts['assemble_wing']
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	faces = f.findAt(((0.193574, 0.015734, 0.0), ))
	leaf = dgm.LeafFromGeometry(faceSeq=faces)
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	faces = f.findAt(((0.218246, 0.053864, 0.003333), ))
	leaf = dgm.LeafFromGeometry(faceSeq=faces)
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	faces = f.findAt(((0.289629, 0.053867, 0.006667), ))
	leaf = dgm.LeafFromGeometry(faceSeq=faces)
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	faces = f.findAt(((0.242298, 0.053867, 0.34), ))
	leaf = dgm.LeafFromGeometry(faceSeq=faces)
	leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	faces = f.findAt(((0.193574, 0.015734, 0.0), ))
	leaf = dgm.LeafFromGeometry(faceSeq=faces)
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	faces = f.findAt(((0.152274, 0.046163, 0.0), ))
	leaf = dgm.LeafFromGeometry(faceSeq=faces)
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	faces = f.findAt(((0.218246, 0.053864, 0.003333), ))
	leaf = dgm.LeafFromGeometry(faceSeq=faces)
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	faces = f.findAt(((0.287576, -0.060397, 0.34), ))
	leaf = dgm.LeafFromGeometry(faceSeq=faces)
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	faces = f.findAt(((0.216898, -0.0604, 1.666667), ))
	leaf = dgm.LeafFromGeometry(faceSeq=faces)
	leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
	a = mdb.models['Model-1'].rootAssembly
	a1 = mdb.models['Model-1'].rootAssembly
	i1 = a1.instances['assemble_wing-1']
	leaf = dgm.LeafFromInstance((i1, ))
	leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
	a1 = mdb.models['Model-1'].rootAssembly
	f1 = a1.instances['Clevice-1'].faces
	faces1 = f1.findAt(((0.196233, -0.053733, 0.0), ), ((0.196233, 0.047198, 0.0), 
		), ((0.196233, -0.03673, 0.0), ), ((0.199567, 0.030194, 0.0), ))
	leaf = dgm.LeafFromGeometry(faceSeq=faces1)
	session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	a1 = mdb.models['Model-1'].rootAssembly
	i1 = a1.instances['Clevice-1']
	leaf = dgm.LeafFromInstance((i1, ))
	a1 = mdb.models['Model-1'].rootAssembly
	f1 = a1.instances['assemble_wing-1'].faces
	faces1 = f1.findAt(((0.152274, 0.046163, 0.0), ), ((0.193574, 0.015734, 0.0), 
		))
	leaf = dgm.LeafFromGeometry(faceSeq=faces1)
	a1 = mdb.models['Model-1'].rootAssembly
	s1 = a1.instances['Clevice-1'].faces
	side1Faces1 = s1.findAt(((0.196233, -0.053733, 0.0), ), ((0.196233, 0.047198, 
		0.0), ), ((0.196233, -0.03673, 0.0), ), ((0.199567, 0.030194, 0.0), ))
	region1=a1.Surface(side1Faces=side1Faces1, name='m_Surf-1')
	a1 = mdb.models['Model-1'].rootAssembly
	s1 = a1.instances['assemble_wing-1'].edges
	side1Edges1 = s1.findAt(((0.2474, 0.053864, 0.0), ), ((0.1929, 0.025198, 0.0), 
		), ((0.210898, -0.0604, 0.0), ))
	region2=a1.Surface(side1Edges=side1Edges1, name='s_Surf-1')
	mdb.models['Model-1'].Tie(name='tie_spar_clevis', master=region1, 
		slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
		tieRotations=ON, thickness=ON)
	leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
    
    # Creating Partitioning Datum for Assembly Mesh Region:
	p = mdb.models['Model-1'].parts['assemble_wing']
	datum = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=selectExclude2-0.0001)
	datumNum = datum.id
	f = p.faces
	pickedFaces = f[:]
	d = p.datums
	p.PartitionFaceByDatumPlane(datumPlane=d[datumNum], faces=pickedFaces)
	# Defining Assembly Mesh Set:
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	faces = f.getByBoundingBox(-10.0,-10.0,selectExclude2-0.001,10.0,10.0,10.0)
	p.Set(faces=faces, name='meshSet2')
	a=mdb.models['Model-1'].rootAssembly
	a.SetByBoolean(name='SCANME', sets=( 
	a.allInstances['Clevice-1'].sets['meshSet1'], ))
	
	
	# a = mdb.models['Model-1'].rootAssembly
	# i1 = a.instances['assemble_wing-1']
	# leaf = dgm.LeafFromInstance((i1, ))
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=12.8745, 
		# farPlane=13.2253, width=0.811596, height=0.375787, viewOffsetX=1.68915, 
		# viewOffsetY=1.19188)
	# session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
	# leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
		# leaf=leaf)
	# session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=9.41946, 
		# farPlane=13.476, width=2.4643, height=1.14103, viewOffsetX=0.758701, 
		# viewOffsetY=0.467289)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=9.33263, 
		# farPlane=14.9394, width=2.44159, height=1.13051, cameraPosition=(1.1959, 
		# 2.75291, -9.45469), cameraUpVector=(-0.928568, 0.0450602, 0.368416), 
		# cameraTarget=(1.50556, 1.16297, 1.87787), viewOffsetX=0.751707, 
		# viewOffsetY=0.462981)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=9.35869, 
		# farPlane=14.9133, width=2.05867, height=0.95321, viewOffsetX=1.37098, 
		# viewOffsetY=0.849731)
	# a = mdb.models['Model-1'].rootAssembly
	# f1 = a.instances['Clevice-1'].faces
	# faces1 = f1.findAt(((0.19791, -0.05038, 0.0), ), ((0.19791, 0.043844, 0.0), ), 
		# ((0.19791, -0.034503, 0.0), ), ((0.20292, 0.027967, 0.0), ))
	# leaf = dgm.LeafFromGeometry(faceSeq=faces1)
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	# a = mdb.models['Model-1'].rootAssembly
	# i1 = a.instances['Clevice-1']
	# leaf = dgm.LeafFromInstance((i1, ))
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=9.57184, 
		# farPlane=14.8351, width=1.45256, height=0.672568, viewOffsetX=1.43037, 
		# viewOffsetY=0.924165)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=9.54877, 
		# farPlane=14.7613, width=1.44906, height=0.670947, cameraPosition=(1.44389, 
		# 3.52237, -9.18391), cameraUpVector=(-0.932028, 0.0446282, 0.359627), 
		# cameraTarget=(1.55842, 1.15995, 2.01683), viewOffsetX=1.42693, 
		# viewOffsetY=0.921937)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=9.4473, 
		# farPlane=14.8628, width=2.85487, height=1.32187, viewOffsetX=1.34149, 
		# viewOffsetY=0.8287)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=9.43969, 
		# farPlane=14.5926, width=2.85257, height=1.3208, cameraPosition=(2.11526, 
		# 6.09928, -7.82898), cameraUpVector=(-0.94813, -0.0121114, 0.317651), 
		# cameraTarget=(1.60961, 1.16078, 2.48636), viewOffsetX=1.34041, 
		# viewOffsetY=0.828032)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=9.63102, 
		# farPlane=14.4013, width=0.79691, height=0.368987, viewOffsetX=2.20245, 
		# viewOffsetY=1.0997)
	# a = mdb.models['Model-1'].rootAssembly
	# f1 = a.instances['assemble_wing-1'].faces
	# faces1 = f1.findAt(((0.193574, 0.015734, 0.0), ), ((0.152274, 0.046163, 0.0), 
		# ), ((0.289629, 0.053867, 0.006667), ))
	# leaf = dgm.LeafFromGeometry(faceSeq=faces1)
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	# a = mdb.models['Model-1'].rootAssembly
	# f1 = a.instances['assemble_wing-1'].faces
	# faces1 = f1.findAt(((0.240276, -0.060397, 0.006667), ))
	# leaf = dgm.LeafFromGeometry(faceSeq=faces1)
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=9.71113, 
		# farPlane=14.3212, width=0.00640759, height=0.00296685, viewOffsetX=2.24651, 
		# viewOffsetY=1.19652)
	# a = mdb.models['Model-1'].rootAssembly
	# s1 = a.instances['Clevice-1'].faces
	# side1Faces1 = s1.findAt(((0.19791, -0.05038, 0.0), ), ((0.19791, 0.043844, 
		# 0.0), ), ((0.19791, -0.034503, 0.0), ), ((0.20292, 0.027967, 0.0), ))
	# region1=a.Surface(side1Faces=side1Faces1, name='m_Surf-14')
	# a = mdb.models['Model-1'].rootAssembly
	# s1 = a.instances['assemble_wing-1'].edges
	# side1Edges1 = s1.findAt(((0.210898, -0.0604, 0.0), ))
	# side2Edges1 = s1.findAt(((0.2474, 0.053864, 0.0), ), ((0.1929, 0.025198, 0.0), 
		# ), ((0.194417, 0.05383, 0.0), ))
	# region2=a.Surface(side1Edges=side1Edges1, side2Edges=side2Edges1, 
		# name='s_Surf-14')
	# mdb.models['Model-1'].Tie(name='tie_spar_clevis', master=region1, 
		# slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
		# tieRotations=ON, thickness=ON)
	# leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
		# leaf=leaf)
	# session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
	
	# # w/o paritions
	# a = mdb.models['Model-1'].rootAssembly
	# i1 = a.instances['assemble_wing-1']
	# leaf = dgm.LeafFromInstance((i1, ))
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=12.8794, 
		# farPlane=13.2204, width=0.76319, height=0.353374, viewOffsetX=1.63063, 
		# viewOffsetY=1.27267)
	# leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
		# leaf=leaf)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=9.20619, 
		# farPlane=13.6893, width=4.29005, height=1.98639, viewOffsetX=1.06645, 
		# viewOffsetY=1.03293)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=10.9567, 
		# farPlane=14.8007, width=5.1058, height=2.3641, cameraPosition=(7.20041, 
		# 9.57419, -3.12273), cameraUpVector=(-0.977314, 0.186726, 0.099948), 
		# cameraTarget=(2.20777, 0.902343, 2.43815), viewOffsetX=1.26923, 
		# viewOffsetY=1.22934)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=11.2007, 
		# farPlane=14.5567, width=1.63834, height=0.758586, viewOffsetX=2.66259, 
		# viewOffsetY=1.42614)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=11.4365, 
		# farPlane=15.8038, width=1.67282, height=0.77455, cameraPosition=(5.89812, 
		# 8.21307, -7.16991), cameraUpVector=(-0.866757, 0.294031, 0.402839), 
		# cameraTarget=(2.85513, 0.89877, 1.094), viewOffsetX=2.71862, 
		# viewOffsetY=1.45615)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=11.4334, 
		# farPlane=15.807, width=1.57202, height=0.727879, viewOffsetX=2.59098, 
		# viewOffsetY=1.42297)
	# a = mdb.models['Model-1'].rootAssembly
	# f1 = a.instances['Clevice-1'].faces
	# faces1 = f1.findAt(((0.20292, -0.023635, 0.0), ), ((0.20292, 0.0171, 0.0), ))
	# leaf = dgm.LeafFromGeometry(faceSeq=faces1)
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	# a = mdb.models['Model-1'].rootAssembly
	# i1 = a.instances['Clevice-1']
	# leaf = dgm.LeafFromInstance((i1, ))
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=11.5196, 
		# farPlane=15.8191, width=1.69187, height=0.783374, viewOffsetX=2.59112, 
		# viewOffsetY=1.50127)
	# a = mdb.models['Model-1'].rootAssembly
	# f1 = a.instances['assemble_wing-1'].faces
	# faces1 = f1.findAt(((0.193574, 0.015734, 0.0), ), ((0.152274, 0.046163, 0.0), 
		# ))
	# leaf = dgm.LeafFromGeometry(faceSeq=faces1)
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=11.6417, 
		# farPlane=15.697, width=0.549677, height=0.254512, viewOffsetX=2.50997, 
		# viewOffsetY=1.44986)
	# a = mdb.models['Model-1'].rootAssembly
	# f1 = a.instances['assemble_wing-1'].faces
	# faces1 = f1.findAt(((0.240276, -0.060397, 0.006667), ), ((0.289629, 0.053867, 
		# 0.006667), ))
	# leaf = dgm.LeafFromGeometry(faceSeq=faces1)
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.remove(leaf=leaf)
	# session.viewports['Viewport: 1'].view.setValues(nearPlane=11.6547, 
		# farPlane=15.6841, width=0.40551, height=0.18776, viewOffsetX=2.51467, 
		# viewOffsetY=1.45548)
	# a = mdb.models['Model-1'].rootAssembly
	# s1 = a.instances['Clevice-1'].faces
	# side1Faces1 = s1.findAt(((0.20292, -0.023635, 0.0), ), ((0.20292, 0.0171, 0.0), 
		# ))
	# region1=a.Surface(side1Faces=side1Faces1, name='m_Surf-14')
	# a = mdb.models['Model-1'].rootAssembly
	# s1 = a.instances['assemble_wing-1'].edges
	# side1Edges1 = s1.findAt(((0.210898, -0.0604, 0.0), ))
	# side2Edges1 = s1.findAt(((0.2474, 0.053864, 0.0), ), ((0.1929, 0.025198, 0.0), 
		# ))
	# region2=a.Surface(side1Edges=side1Edges1, side2Edges=side2Edges1, 
		# name='s_Surf-14')
	# mdb.models['Model-1'].Tie(name='tie_spar_clevis', master=region1, 
		# slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
		# tieRotations=ON, thickness=ON)
	# leaf = dgm.Leaf(leafType=DEFAULT_MODEL)
	# session.viewports['Viewport: 1'].assemblyDisplay.displayGroup.replace(
		# leaf=leaf)
	

	#Mesh Parts
	print 'Meshing the Part'
	# wing assemble
	p = mdb.models['Model-1'].parts['assemble_wing']
	f = p.faces
	#pickedRegions = f[:]
	pickedRegions = f.findAt(((0.388651, -0.061133, 0.006667), ), ((0.388651, 
		-0.061133, 0.67), ), ((0.388651, -0.061133, 1.006667), ), ((0.388651, 
		-0.061133, 1.67), ), ((0.388651, -0.061133, 2.006667), ), ((0.388651, 
		-0.061133, 2.67), ), ((0.388651, -0.061133, 3.006667), ), ((0.388651, 
		-0.061133, 3.67), ), ((0.388651, -0.061133, 4.006667), ), ((0.388651, 
		-0.061133, 4.67), ), ((0.481696, -0.060507, 4.003333), ), ((0.481696, 
		-0.060507, 3.003333), ), ((0.481696, -0.060507, 2.003333), ), ((0.481696, 
		-0.060507, 1.003333), ), ((0.481696, -0.060507, 0.67), ), ((0.536253, 
		-0.059323, 0.006667), ), ((0.536253, -0.059323, 0.67), ), ((0.536253, 
		-0.059323, 1.006667), ), ((0.536253, -0.059323, 1.67), ), ((0.536253, 
		-0.059323, 2.006667), ), ((0.536253, -0.059323, 2.67), ), ((0.536253, 
		-0.059323, 3.006667), ), ((0.536253, -0.059323, 3.67), ), ((0.536253, 
		-0.059323, 4.006667), ), ((0.536253, -0.059323, 4.67), ), ((0.630622, 
		-0.052068, 4.003333), ), ((0.630622, -0.052068, 3.003333), ), ((0.630622, 
		-0.052068, 2.003333), ), ((0.630622, -0.052068, 1.003333), ), ((0.630622, 
		-0.052068, 0.67), ), ((0.685453, -0.049569, 0.006667), ), ((0.685453, 
		-0.049569, 0.67), ), ((0.685453, -0.049569, 1.006667), ), ((0.685453, 
		-0.049569, 1.67), ), ((0.685453, -0.049569, 2.006667), ), ((0.685453, 
		-0.049569, 2.67), ), ((0.685453, -0.049569, 3.006667), ), ((0.685453, 
		-0.049569, 3.67), ), ((0.685453, -0.049569, 4.006667), ), ((0.685453, 
		-0.049569, 4.67), ), ((0.780267, -0.040032, 4.003333), ), ((0.780267, 
		-0.040032, 3.003333), ), ((0.780267, -0.040032, 2.003333), ), ((0.780267, 
		-0.040032, 1.003333), ), ((0.780267, -0.040032, 0.67), ), ((0.83441, 
		-0.036501, 0.006667), ), ((0.83441, -0.036501, 0.67), ), ((0.83441, 
		-0.036501, 1.006667), ), ((0.83441, -0.036501, 1.67), ), ((0.83441, 
		-0.036501, 2.006667), ), ((0.83441, -0.036501, 2.67), ), ((0.83441, 
		-0.036501, 3.006667), ), ((0.83441, -0.036501, 3.67), ), ((0.83441, 
		-0.036501, 4.006667), ), ((0.83441, -0.036501, 4.67), ), ((0.927768, 
		-0.025198, 0.67), ), ((0.927768, -0.025198, 1.67), ), ((0.927768, 
		-0.025198, 2.67), ), ((0.927768, -0.025198, 4.003333), ), ((0.927768, 
		-0.025198, 3.67), ), ((0.932756, -0.026534, 4.34), ), ((0.932756, 
		-0.026534, 4.003333), ), ((0.932756, -0.026534, 3.34), ), ((0.932756, 
		-0.026534, 3.003333), ), ((0.932756, -0.026534, 2.34), ), ((0.932756, 
		-0.026534, 2.003333), ), ((0.932756, -0.026534, 1.34), ), ((0.932756, 
		-0.026534, 1.003333), ), ((0.932756, -0.026534, 0.34), ), ((0.932756, 
		-0.026534, 0.003333), ), ((0.336826, 0.051471, 4.003333), ), ((0.336826, 
		0.051471, 3.003333), ), ((0.336826, 0.051471, 2.003333), ), ((0.336826, 
		0.051471, 0.67), ), ((0.336826, 0.051471, 1.003333), ), ((0.340764, 
		0.053652, 0.006667), ), ((0.340764, 0.053652, 0.67), ), ((0.340764, 
		0.053652, 1.006667), ), ((0.340764, 0.053652, 1.67), ), ((0.340764, 
		0.053652, 2.006667), ), ((0.340764, 0.053652, 2.67), ), ((0.340764, 
		0.053652, 3.006667), ), ((0.340764, 0.053652, 3.67), ), ((0.340764, 
		0.053652, 4.006667), ), ((0.340764, 0.053652, 4.67), ), ((0.483469, 
		0.043166, 4.003333), ), ((0.483469, 0.043166, 3.003333), ), ((0.483469, 
		0.043166, 2.003333), ), ((0.483469, 0.043166, 0.67), ), ((0.483469, 
		0.043166, 1.003333), ), ((0.487482, 0.045207, 0.006667), ), ((0.487482, 
		0.045207, 0.67), ), ((0.487482, 0.045207, 1.006667), ), ((0.487482, 
		0.045207, 1.67), ), ((0.487482, 0.045207, 2.006667), ), ((0.487482, 
		0.045207, 2.67), ), ((0.487482, 0.045207, 3.006667), ), ((0.487482, 
		0.045207, 3.67), ), ((0.487482, 0.045207, 4.006667), ), ((0.487482, 
		0.045207, 4.67), ), ((0.631998, 0.029536, 4.003333), ), ((0.631998, 
		0.029536, 3.003333), ), ((0.631998, 0.029536, 2.003333), ), ((0.631998, 
		0.029536, 0.67), ), ((0.631998, 0.029536, 1.003333), ), ((0.636057, 
		0.031482, 0.006667), ), ((0.636057, 0.031482, 0.67), ), ((0.636057, 
		0.031482, 1.006667), ), ((0.636057, 0.031482, 1.67), ), ((0.636057, 
		0.031482, 2.006667), ), ((0.636057, 0.031482, 2.67), ), ((0.636057, 
		0.031482, 3.006667), ), ((0.636057, 0.031482, 3.67), ), ((0.636057, 
		0.031482, 4.006667), ), ((0.636057, 0.031482, 4.67), ), ((0.781122, 
		0.012286, 4.003333), ), ((0.781122, 0.012286, 3.003333), ), ((0.781122, 
		0.012286, 2.003333), ), ((0.781122, 0.012286, 0.67), ), ((0.781122, 
		0.012286, 1.003333), ), ((0.785219, 0.014151, 0.006667), ), ((0.785219, 
		0.014151, 0.67), ), ((0.785219, 0.014151, 1.006667), ), ((0.785219, 
		0.014151, 1.67), ), ((0.785219, 0.014151, 2.006667), ), ((0.785219, 
		0.014151, 2.67), ), ((0.785219, 0.014151, 3.006667), ), ((0.785219, 
		0.014151, 3.67), ), ((0.785219, 0.014151, 4.006667), ), ((0.785219, 
		0.014151, 4.67), ), ((0.928008, -0.007686, 4.003333), ), ((0.928008, 
		-0.007686, 3.003333), ), ((0.928008, -0.007686, 2.003333), ), ((0.928008, 
		-0.007686, 0.67), ), ((0.928008, -0.007686, 1.003333), ), ((0.932134, 
		-0.005886, 0.006667), ), ((0.932134, -0.005886, 0.67), ), ((0.932134, 
		-0.005886, 1.006667), ), ((0.932134, -0.005886, 1.67), ), ((0.932134, 
		-0.005886, 2.006667), ), ((0.932134, -0.005886, 2.67), ), ((0.932134, 
		-0.005886, 3.006667), ), ((0.932134, -0.005886, 3.67), ), ((0.932134, 
		-0.005886, 4.006667), ), ((0.932134, -0.005886, 4.67), ), ((0.218246, 
		0.053864, 1.003333), ), ((0.218246, 0.053864, 2.003333), ), ((0.218246, 
		0.053864, 3.003333), ), ((0.218246, 0.053864, 4.003333), ), ((0.218246, 
		0.053864, 0.003333), ), ((0.1929, -0.022351, 1.67), ), ((0.1929, 0.015689, 
		0.003333), ), ((0.1929, -0.060393, 0.003333), ), ((0.1929, -0.022351, 
		2.67), ), ((0.1929, 0.053773, 0.003333), ), ((0.1929, -0.022351, 4.67), ), 
		((0.194248, 0.053864, 0.003333), ), ((0.1929, -0.022351, 3.67), ), ((
		0.216898, -0.0604, 1.666667), ), ((0.192925, -0.060392, 0.006667), ), ((
		0.158159, 0.051232, 0.006667), ), ((0.194248, 0.053819, 0.006667), ), ((
		0.192925, -0.060392, 0.67), ), ((0.158159, 0.051232, 0.67), ), ((0.194248, 
		0.053819, 0.67), ), ((0.192925, -0.060392, 1.006667), ), ((0.158159, 
		0.051232, 1.006667), ), ((0.194248, 0.053819, 1.006667), ), ((0.192925, 
		-0.060392, 1.67), ), ((0.158159, 0.051232, 1.67), ), ((0.194248, 0.053819, 
		1.67), ), ((0.192925, -0.060392, 2.006667), ), ((0.158159, 0.051232, 
		2.006667), ), ((0.194248, 0.053819, 2.006667), ), ((0.192925, -0.060392, 
		2.67), ), ((0.158159, 0.051232, 2.67), ), ((0.194248, 0.053819, 2.67), ), (
		(0.192925, -0.060392, 3.006667), ), ((0.158159, 0.051232, 3.006667), ), ((
		0.194248, 0.053819, 3.006667), ), ((0.192925, -0.060392, 3.67), ), ((
		0.158159, 0.051232, 3.67), ), ((0.194248, 0.053819, 3.67), ), ((0.192925, 
		-0.060392, 4.006667), ), ((0.158159, 0.051232, 4.006667), ), ((0.194248, 
		0.053819, 4.006667), ), ((0.192925, -0.060392, 4.67), ), ((0.158159, 
		0.051232, 4.67), ), ((0.194248, 0.053819, 4.67), ), ((0.242298, 0.053867, 
		0.34), ), ((0.194937, 0.053865, 0.003333), ), ((0.389673, 0.050881, 0.34), 
		), ((0.240276, -0.060397, 0.006667), ), ((0.53701, 0.040662, 0.34), ), ((
		0.339675, -0.060463, 0.003333), ), ((0.685781, 0.02573, 0.34), ), ((
		0.486596, -0.062137, 0.003333), ), ((0.834194, 0.007492, 0.34), ), ((
		0.635558, -0.053582, 0.003333), ), ((0.955968, -0.009513, 0.34), ), ((
		0.785233, -0.041447, 0.003333), ), ((0.999855, -0.017032, 0.34), ), ((
		0.956702, -0.023741, 0.67), ), ((0.242298, 0.053867, 1.34), ), ((0.242298, 
		0.053867, 1.003333), ), ((0.194937, 0.053865, 1.34), ), ((0.194937, 
		0.053865, 0.34), ), ((0.389673, 0.050881, 1.34), ), ((0.389673, 0.050881, 
		1.003333), ), ((0.287576, -0.060397, 1.34), ), ((0.287576, -0.060397, 
		0.34), ), ((0.53701, 0.040662, 1.34), ), ((0.53701, 0.040662, 1.003333), ), 
		((0.339675, -0.060463, 1.34), ), ((0.339675, -0.060463, 0.34), ), ((
		0.685781, 0.02573, 1.34), ), ((0.685781, 0.02573, 1.003333), ), ((0.486596, 
		-0.062137, 1.34), ), ((0.486596, -0.062137, 0.34), ), ((0.834194, 0.007492, 
		1.34), ), ((0.834194, 0.007492, 1.003333), ), ((0.635558, -0.053582, 1.34), 
		), ((0.635558, -0.053582, 0.34), ), ((0.955968, -0.009513, 1.34), ), ((
		0.955968, -0.009513, 1.003333), ), ((0.785233, -0.041447, 1.34), ), ((
		0.785233, -0.041447, 0.34), ), ((0.999855, -0.017032, 1.34), ), ((0.999855, 
		-0.017032, 1.003333), ), ((0.956702, -0.023741, 1.67), ), ((0.956702, 
		-0.023741, 1.006667), ), ((0.242298, 0.053867, 2.34), ), ((0.242298, 
		0.053867, 2.003333), ), ((0.194937, 0.053865, 2.34), ), ((0.194937, 
		0.053865, 1.003333), ), ((0.389673, 0.050881, 2.34), ), ((0.389673, 
		0.050881, 2.003333), ), ((0.287576, -0.060397, 2.34), ), ((0.287576, 
		-0.060397, 1.003333), ), ((0.53701, 0.040662, 2.34), ), ((0.53701, 
		0.040662, 2.003333), ), ((0.339675, -0.060463, 2.34), ), ((0.339675, 
		-0.060463, 1.003333), ), ((0.685781, 0.02573, 2.34), ), ((0.685781, 
		0.02573, 2.003333), ), ((0.486596, -0.062137, 2.34), ), ((0.486596, 
		-0.062137, 1.003333), ), ((0.834194, 0.007492, 2.34), ), ((0.834194, 
		0.007492, 2.003333), ), ((0.635558, -0.053582, 2.34), ), ((0.635558, 
		-0.053582, 1.003333), ), ((0.955968, -0.009513, 2.34), ), ((0.955968, 
		-0.009513, 2.003333), ), ((0.785233, -0.041447, 2.34), ), ((0.785233, 
		-0.041447, 1.003333), ))+f.findAt(((0.999855, -0.017032, 2.34), ), ((
		0.999855, -0.017032, 2.003333), ), ((0.956702, -0.023741, 2.67), ), ((
		0.956702, -0.023741, 2.006667), ), ((0.242298, 0.053867, 3.34), ), ((
		0.242298, 0.053867, 3.003333), ), ((0.194937, 0.053865, 3.34), ), ((
		0.194937, 0.053865, 2.003333), ), ((0.389673, 0.050881, 3.34), ), ((
		0.389673, 0.050881, 3.003333), ), ((0.287576, -0.060397, 3.34), ), ((
		0.287576, -0.060397, 2.003333), ), ((0.53701, 0.040662, 3.34), ), ((
		0.53701, 0.040662, 3.003333), ), ((0.339675, -0.060463, 3.34), ), ((
		0.339675, -0.060463, 2.003333), ), ((0.685781, 0.02573, 3.34), ), ((
		0.685781, 0.02573, 3.003333), ), ((0.486596, -0.062137, 3.34), ), ((
		0.486596, -0.062137, 2.003333), ), ((0.834194, 0.007492, 3.34), ), ((
		0.834194, 0.007492, 3.003333), ), ((0.635558, -0.053582, 3.34), ), ((
		0.635558, -0.053582, 2.003333), ), ((0.955968, -0.009513, 3.34), ), ((
		0.955968, -0.009513, 3.003333), ), ((0.785233, -0.041447, 3.34), ), ((
		0.785233, -0.041447, 2.003333), ), ((0.999855, -0.017032, 3.34), ), ((
		0.999855, -0.017032, 3.003333), ), ((0.956702, -0.023741, 3.67), ), ((
		0.956702, -0.023741, 3.006667), ), ((0.242298, 0.053867, 4.34), ), ((
		0.242298, 0.053867, 4.003333), ), ((0.194937, 0.053865, 4.34), ), ((
		0.194937, 0.053865, 3.003333), ), ((0.389673, 0.050881, 4.34), ), ((
		0.389673, 0.050881, 4.003333), ), ((0.287576, -0.060397, 4.34), ), ((
		0.287576, -0.060397, 3.003333), ), ((0.53701, 0.040662, 4.34), ), ((
		0.53701, 0.040662, 4.003333), ), ((0.339675, -0.060463, 4.34), ), ((
		0.339675, -0.060463, 3.003333), ), ((0.685781, 0.02573, 4.34), ), ((
		0.685781, 0.02573, 4.003333), ), ((0.486596, -0.062137, 4.34), ), ((
		0.486596, -0.062137, 3.003333), ), ((0.834194, 0.007492, 4.34), ), ((
		0.834194, 0.007492, 4.003333), ), ((0.635558, -0.053582, 4.34), ), ((
		0.635558, -0.053582, 3.003333), ), ((0.955968, -0.009513, 4.34), ), ((
		0.955968, -0.009513, 4.003333), ), ((0.785233, -0.041447, 4.34), ), ((
		0.785233, -0.041447, 3.003333), ), ((0.999855, -0.017032, 4.34), ), ((
		0.999855, -0.017032, 4.003333), ), ((0.956702, -0.023741, 4.67), ), ((
		0.956702, -0.023741, 4.006667), ), ((0.978264, -0.021227, 5.003333), ), ((
		0.881199, -0.031795, 5.003333), ), ((0.732955, -0.045748, 5.003333), ), ((
		0.583513, -0.056645, 5.003333), ), ((0.435228, -0.061771, 5.003333), ), ((
		0.287576, -0.060397, 5.006667), ), ((0.022632, -0.025555, 5.006667), ), ((
		0.242298, 0.053867, 5.006667), ), ((0.338863, 0.05376, 5.003333), ), ((
		0.485585, 0.045381, 5.003333), ), ((0.634165, 0.031701, 5.003333), ), ((
		0.783332, 0.014408, 5.003333), ), ((0.930251, -0.0056, 5.003333), ), ((
		0.999855, -0.017032, 5.006667), ), ((0.955968, -0.009513, 0.003333), ), ((
		0.834194, 0.007492, 0.003333), ), ((0.685781, 0.02573, 0.003333), ), ((
		0.53701, 0.040662, 0.003333), ), ((0.389673, 0.050881, 0.003333), ), ((
		0.289629, 0.053867, 0.006667), ), ((0.194937, 0.053865, 4.003333), ), ((
		0.287576, -0.060397, 4.003333), ), ((0.339675, -0.060463, 4.003333), ), ((
		0.486596, -0.062137, 4.003333), ), ((0.635558, -0.053582, 4.003333), ), ((
		0.785233, -0.041447, 4.003333), ), ((0.956702, -0.023741, 0.006667), ), ((
		0.99984, -0.017872, 0.006667), ))
	p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
	p = mdb.models['Model-1'].parts['assemble_wing']
	p.seedPart(size=assemblySeed, deviationFactor=0.1, minSizeFactor=0.1)
	p = mdb.models['Model-1'].parts['assemble_wing']
	p.generateMesh()

	
	# Clevis
	p = mdb.models['Model-1'].parts['Clevice']
	p = mdb.models['Model-1'].parts['Clevice']
	p.seedPart(size=clevisSeed, deviationFactor=0.1, minSizeFactor=0.1)
	p = mdb.models['Model-1'].parts['Clevice']
	p.generateMesh()
	
	

	#####################################
	### Creation/Execution of the Job ###
	#####################################
	print 'Creating/Running Job'

	ModelName='Model-1'

	mdb.Job(name=ModelName, model=ModelName, description='', type=ANALYSIS, 
		atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
		memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
		explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
		modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
		scratch='', multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)

	job=mdb.jobs[ModelName]

	# delete lock file, which for some reason tends to hang around, if it exists
	if os.access('%s.lck'%ModelName,os.F_OK):
		os.remove('%s.lck'%ModelName)
		
	# Run the job, then process the results.		
	job.submit()
	job.waitForCompletion()
	print 'Completed job'
	
	# Finds Eigen Value and Checks for Buckling
	eig = findEigenValue('Model-1','Buckle-Step')
	crit = 1.0 # need to figure out what this is
	print eig
	if (eig <= crit):
		buckle = 1
		print 'Buckled'
	else:
		buckle = 0
		print 'Did Not Buckle'
	print 'Buckle or not: %d (1: buckle; 0: no buckle)' %(buckle)
	
	pr=mdb.models[ModelName].rootAssembly.getMassProperties()
	mass=pr['mass']
	print 'Total mass: %f kg' %(mass)
	
	print 'Made it in'
	odbName = ModelName+'.odb'
	print odbName
	odb = visualization.openOdb(odbName)
	lastFrame = odb.steps['Load-Step'].frames[-1]
	print 'odb open'

	## The following is left for use in later probems/projects
	print 'Scanning the PART for maximum VM STRESS'
	elsetName='SCANME'
	elset = elemset = None
	region = "over the entire model"
	assembly = odb.rootAssembly

	#Check to see if the element set exists
	#in the assembly

	if elsetName:
		try:
			elemset = assembly.elementSets[elsetName]
			region = " in the element set : " + elsetName;
		except KeyError:
			print 'An assembly level elset named %s does' \
				  'not exist in the output database %s' \
				  % (elsetName, odbName)			
			
	maxMises = -0.1
	maxVMElem = 0
	maxStep = "_None_"
	maxFrame = -1
	Stress = 'S'
	Disp = 'U'
	isStressPresent = 0
	for step in odb.steps.values():
		print 'Processing Step:', step.name
		for frame in step.frames:
			allFields = frame.fieldOutputs
			if (allFields.has_key(Stress)):
				isStressPresent = 1
				stressSet = allFields[Stress]
				if elemset:
					stressSet = stressSet.getSubset(
						region=elemset)		 
				for stressValue in stressSet.values:				
					if (stressValue.mises > maxMises):
						maxMises = stressValue.mises
						maxVMElem = stressValue.elementLabel
						maxStep = step.name
						maxFrame = frame.incrementNumber
	if(isStressPresent):
		print 'Maximum von Mises stress %s is %f in element %d'%(
			region, maxMises, maxVMElem)
		print 'Location: frame # %d	 step:	%s '%(maxFrame,maxStep)
	else:
		print 'Stress output is not available in' \
			  'the output database : %s\n' %(odb.name)
	
	
	if (maxMises >= YieldTensile):
		print 'Material Yielded'
		Yield = 1
	else:
		print 'Material did not Yield'
		Yield = 0
	
	if (buckle == 0 and	 Yield == 0):
		print 'Design Success'
		Design = 0
	else:
		print 'Design Failure'
		Design = 1
	
	print 'Retrieving ALL displacements at TIPNODE'
	TIP_SET_RIB = 'RIB-'+str(numRibs-1)
	pTip = odb.rootAssembly.nodeSets['TIPNODE']
	
	# Retrieve Y-displacements at the splines/connectors
	print 'Retrieving ALL final displacements at ALL points'
	dispField = lastFrame.fieldOutputs['U']

	dFieldpTip = dispField.getSubset(region=pTip)
	
	print 'Retrieving only U3 at TIPNODE'
	#Note, U1=data[0], U2=data[1], U3=data[2]
	disppTipNode = dFieldpTip.values[0].data[1]
	
	# Extract displacements
	# Selecting the node(s) to be queried
	pTip = odb.rootAssembly.nodeSets['RIB-OML-0']
			
	# Retrieve Y-displacements at the splines/connectors
	print 'Retrieving ALL final displacements at ALL points'
	dispField = lastFrame.fieldOutputs['U']

	print 'Retrieving ALL displacements at NODE'
	dFieldpTip = dispField.getSubset(region=pTip)
		
	print 'Retrieving only U2 at NODE on the rib'
	#Note, U1=data[0], U2=data[1], U3=data[2]
	disppTipX = [[0]*len(dFieldpTip.values)]*numRibs
	disppTipY = [[0]*len(dFieldpTip.values)]*numRibs
	for j in range(numRibs):
		for i in range(len(dFieldpTip.values)):
			disppTipX[j][i] = dFieldpTip.values[i].data[0]
			disppTipY[j][i] = dFieldpTip.values[i].data[1]
			# print disppTipX[i]
			# print disppTipY[i]

	DataFile = open(filenameDisppTipX, 'a')
	for j in range(numRibs):
		for i in range(len(dFieldpTip.values)):
			DataFile.write('%f\n' %(disppTipX[j][i]))
	DataFile.close()
	DataFile = open(filenameDisppTipY, 'a')
	for j in range(numRibs):
		for i in range(len(dFieldpTip.values)):
			DataFile.write('%f\n' %(disppTipY[j][i]))
	DataFile.close()
	odb.close()

	toc = time.time()
	elapsedTime = toc - tic
	print 'The elapsed time is %f sec' %(elapsedTime)
	
	DataFile = open(filename, 'a')
	# DataFile.write('%d %4f %4f %4f %4f %4f %s' %(numRibs, ribThick, SparThick, skinThick, stringThick, csf, MaterialSelection))
	DataFile.write('%10f %10f %10f %10f %10f %10f %10f\n' %(maxMises, mass, disppTipNode, eig, buckle, Yield, Design))
	# DataFile.write('elapsed time = %10f seconds\n' % elapsedTime)
	DataFile.close()
	
	
			

print 'DONE!!'

