"""
Description:
This script solves the ABAQUS portion of HW4. It takes in different
values for mesh size and element types and writes to a file in order to
view the largest tip displacement for each mesh size and element type.

Script by:
Gregory Wilson
Dept. of Aerospace Engineering
Texas A&M University
October 4, 2017 (My birthday)
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
from math import atan, sin, cos, tan
from DOEmethods import LHS
from Post_P_Script_v2 import getResults

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

filename = 'abaqus_results.txt'
### Write data file column headings
DataFile = open(filename,'w')
DataFile.close()

###############################
### Generation of FEA Model ###
###############################

### Note: If you create a loop, START it here 
### (i.e., the contents of the loop below are intended)
# seedSize = [0.25, 0.2, 0.15, 0.1, 0.075]
sizes = 0.05

# meshTypesArray = ['CPE3','CPE6','CPE4','CPE4I','CPE8','CPE8R']
tic = time.time()

h0 = 0.5
L = 1
# # xDesign = [hL,t,l1,l2,l3,y_bottom]
# lb = [0.05,2.5,0.225,0.625,0.85,0.04]
# ub = [0.5,5,0.3,0.675,0.875,0.055]
# nV = len(lb)
# delta_x = [0]*nV
# xm = [0]*nV
# for i in range(nV):
	# delta_x[i] = ub[i]-lb[i]
	# xm[i] = (ub[i]+lb[i])/2

nV = 5 # number of design variables
DV_file = 'DesignVariables.txt'
### Write data file column headings

xDesign = [[0]*nV]
i = 0
file_in = open(DV_file, 'r')
for y in file_in.read().split('\n'):
	if i<nV:
		xDesign[0][i] = (float(y))
		i = i+1
		

# file_in.close()
# t_stiff, h_stiff, w_stiff, n_stiff, n_lam	 = A

# # LHS
# nS = 50
# LHS_samples = LHS(nV,nS)
# xDesign = [[0]*nV]*nS
# for i in range(nS):
	# xDesign_tmp = [0]*nV
	# for j in range(nV):
		# xDesign_tmp[j] = LHS_samples[i][j]*delta_x[j]+lb[j]
	# xDesign[i] = xDesign_tmp

## Full Factorial
# nS = 2**nV
# xDesign = [[0]*nV]*nS
# k = 0
# for i1 in range(2):
	# xDesign_tmp = [0]*nV
	# if i1==0:
		# hL = lb[0]
	# else:
		# hL = ub[0]
	# for i2 in range(2):
		# if i2==0:
			# t = lb[1]
		# else:
			# t = ub[1]
		# for i3 in range(2):
			# if i3==0:
				# l1 = lb[2]
			# else:
				# l1 = ub[2]
			# for i4 in range(2):
				# if i4==0:
					# l2 = lb[3]
				# else:
					# l2 = ub[3]
				# for i5 in range(2):
					# if i5==0:
						# l3 = lb[4]
					# else:
						# l3 = ub[4]
					# for i6 in range(2):
						# if i6==0:
							# y_bot = lb[5]
						# else:
							# y_bot = ub[5]
						# xDesign[k] = [hL,t,l1,l2,l3,y_bot]
						# k = k+1
# print xDesign

## Taguchi
# nS = 18
# Orth_array_file = 'Orth_array_L18.txt'
# ### Write data file column headings
# DataFile = open(Orth_array_file,'r')
# A = [[0]*nV]*nS
# contents = []
# contents = DataFile.read()
# for i in range(nS):
	# A_tmp = [0]*nV
	# for j in range(nV):
		# A_tmp[j] = int(contents[i*12+j*2])
	# A[i] = A_tmp
	# # print(A[i])
# DataFile.close()

# xArray = [lb,xm,ub]
# xDesign = [[0]*nV]*nS
# for i in range(nS):
	# xDesign_tmp = [0]*nV
	# for j in range(nV):
		# index_info = A[i][j]-1
		# xDesign_tmp[j] = xArray[index_info][j]
	# xDesign[i] = xDesign_tmp
	
# print xDesign

## Optimal design
nS = 1
# xDesign = [[0.25000, 2.500000, 0.25, 0.75, 0.1, 0.05]]
for i in range(nS):
	print i
	# DVs
	hL = xDesign[i][0]
	t = xDesign[i][1]
	l1 = 0.25
	l2 = 0.75
	upper_line_slope = (hL-h0)/L
	center_line_slope = upper_line_slope/2
	h_L_over4 = upper_line_slope*l1+h0
	h_L_3over4 = upper_line_slope*l2+h0
	y_l1 = center_line_slope*l1+h0/2
	y_l2 = center_line_slope*l2+h0/2

	# radii
	r1 = xDesign[i][2]*h_L_over4/2
	r2 = xDesign[i][3]*h_L_3over4/2
	# r3 = r1*0.2;
	y1 = y_l1+r1
	y2 = y_l2+r2

	# first point
	# x1 = [l1,l2,l3];
	# y1 = [y_l1,y_l2,y_l3];

	# # second point
	# x2 = [l1,l2,l3];
	# y2 = [y1[0]+r1,y1[1]+r2,y1[2]+r3];

	# # areas of circles
	A1 = pi*r1**2
	A2 = pi*r2**2
	# A3 = pi*r3**2

	# Calculate volume
	volume = t*( (h0+hL)*L/2 - A1 - A2)

	shear_stress = xDesign[i][4]*1e6/(t*hL)

	### Scripting the entire model allows its entire
	### contents to be packaged into this single file. 
	  
	# Sketch Geometry and Create Parts
	print 'Sketching/Creating the part'
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=10.0)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=STANDALONE)
	s.Spot(point=(0.0, 0.0))
	s.Spot(point=(L, 0.0))
	s.Spot(point=(0.0, h0))
	s.Spot(point=(L, hL))
	s.Line(point1=(0.0, h0), point2=(0.0, 0.0))
	s.VerticalConstraint(entity=g[2], addUndoState=False)
	s.Line(point1=(0.0, 0.0), point2=(L, 0.0))
	s.HorizontalConstraint(entity=g[3], addUndoState=False)
	s.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
	s.Line(point1=(L, 0.0), point2=(L, hL))
	s.VerticalConstraint(entity=g[4], addUndoState=False)
	s.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
	s.Line(point1=(L, hL), point2=(0.0, h0))

	s.CircleByCenterPerimeter(center=(l1, y_l1), point1=(l1, y1))
	s.CircleByCenterPerimeter(center=(l2, y_l2), point1=(l2, y2))
	# s.CircleByCenterPerimeter(center=(l3, y_l3), point1=(l3, y2[2]))
	
	p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=TWO_D_PLANAR, 
		type=DEFORMABLE_BODY)
	p = mdb.models['Model-1'].parts['Part-1']
	p.BaseShell(sketch=s)
	s.unsetPrimaryObject()
	p = mdb.models['Model-1'].parts['Part-1']
	session.viewports['Viewport: 1'].setValues(displayedObject=p)
	del mdb.models['Model-1'].sketches['__profile__']
	
	

	#Defining the face partitions
	print 'Partitioning part'
	# p = mdb.models['Model-1'].parts['Part-1']
	# f = p.faces
	# pickedFaces = f.findAt(((0.909958, 0.268907, 0.0), ))
	# v, e, d = p.vertices, p.edges, p.datums
	# p.PartitionFaceByShortestPath(faces=pickedFaces, point1=p.InterestingPoint(
		# edge=e.findAt(coordinates=(0.0, 0.375, 0.0)), rule=MIDDLE), 
		# point2=p.InterestingPoint(edge=e.findAt(coordinates=(L, 0.001, 0.0)), 
		# rule=MIDDLE))
	# p = mdb.models['Model-1'].parts['Part-1']
	# p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=l1)
	# p = mdb.models['Model-1'].parts['Part-1']
	# p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=l2)
	# p = mdb.models['Model-1'].parts['Part-1']
	# p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=l3)
	# p = mdb.models['Model-1'].parts['Part-1']
	# f = p.faces
	# pickedFaces = f.findAt(((0.959319, 0.152892, 0.0), ), ((0.951273, 0.130617, 
		# 0.0), ))
	# d1 = p.datums
	# p.PartitionFaceByDatumPlane(datumPlane=d1[3], faces=pickedFaces)
	# p = mdb.models['Model-1'].parts['Part-1']
	# f = p.faces
	# pickedFaces = f.findAt(((0.283667, 0.404456, 0.0), ), ((0.284212, 0.056379, 
		# 0.0), ), ((0.300365, 0.407249, 0.0), ), ((0.317096, 0.05645, 0.0), ))
	# d = p.datums
	# p.PartitionFaceByDatumPlane(datumPlane=d[4], faces=pickedFaces)
	# p = mdb.models['Model-1'].parts['Part-1']
	# f = p.faces
	# pickedFaces = f.findAt(((0.587155, 0.336014, 0.0), ), ((0.58538, 0.057017, 
		# 0.0), ), ((0.283667, 0.404456, 0.0), ), ((0.284212, 0.056379, 0.0), ), ((
		# 0.600878, 0.33854, 0.0), ), ((0.613216, 0.055554, 0.0), ))
	# d1 = p.datums
	# p.PartitionFaceByDatumPlane(datumPlane=d1[5], faces=pickedFaces)

	# Create Material
	print 'Creating the Materials'
	mdb.models['Model-1'].Material(name='Aluminum')
	mdb.models['Model-1'].materials['Aluminum'].Elastic(table=((70e9, 0.3), 
		))
	
	#Create/Assign Section
	print 'Creating the Sections'
	print 'Assigning the Sections'
	mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', 
	material='Aluminum', thickness=t)
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	faces = f.findAt(((0.01, 0.01, 0.0), ))
	region = p.Set(faces=faces, name='Set-1')
	p = mdb.models['Model-1'].parts['Part-1']
	p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
		offsetType=MIDDLE_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)

	#Assemble Parts
	print 'Placing Parts in Space'
	#Create Instances here
	a = mdb.models['Model-1'].rootAssembly
	session.viewports['Viewport: 1'].setValues(displayedObject=a)
	session.viewports['Viewport: 1'].assemblyDisplay.setValues(
		optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)
	p = mdb.models['Model-1'].parts['Part-1']
	a.Instance(name='Part-1-1', part=p, dependent=ON)

	#Define Steps
	print 'Defining Steps'
	session.viewports['Viewport: 1'].assemblyDisplay.setValues(
		adaptiveMeshConstraints=ON)
	mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', 
		initialInc=0.1)
	session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
	session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
		predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)

	#Create Loads
	print 'Defining Loads'
	v11 = a.instances['Part-1-1'].vertices
	v21 = a.instances['Part-1-1'].vertices
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Part-1-1'].edges
	side1Edges1 = s1.findAt(((1.0, 0.01, 0.0), ))
	region = a.Surface(side1Edges=side1Edges1, name='Surf-1')
	mdb.models['Model-1'].SurfaceTraction(name='Load-1', createStepName='Step-1', 
		region=region, magnitude=-shear_stress, directionVector=(v11.findAt(
		coordinates=(L, 0.0, 0.0)), v21.findAt(coordinates=(L, hL, 0.0))), 
		distributionType=UNIFORM, field='', localCsys=None)
	
	#Define BCs
	print 'Defining all BCs'
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.instances['Part-1-1'].edges
	edges1 = e1.findAt(((0.0, 0.375, 0.0), ))
	region = a.Set(edges=edges1, name='Set-1')
	mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Step-1', 
		region=region, localCsys=None)
	a = mdb.models['Model-1'].rootAssembly
	v1 = a.instances['Part-1-1'].vertices

	#Define Sets
	print 'Defining Sets'
	# Note: Create "TIPNODE" here
	a = mdb.models['Model-1'].rootAssembly
	v1 = a.instances['Part-1-1'].vertices
	verts1 = v1.findAt(((1.0, 0.0, 0.0), ))
	a.Set(vertices=verts1, name='TIPNODE')
	#: The set 'TIPNODE' has been created (1 vertex).
	# Create it from the Assembly Module (not Part)	

	#Mesh Parts
	print 'Meshing the Part'
	p = mdb.models['Model-1'].parts['Part-1']
	p.seedPart(size=sizes, deviationFactor=0.1, minSizeFactor=0.1)
	elemType1 = mesh.ElemType(elemCode=CPS8R, elemLibrary=STANDARD)
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	faces = f.findAt(((0.01, 0.01, 0.0), ))
	pickedRegions =(faces, )
	p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))
	p = mdb.models['Model-1'].parts['Part-1']
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
	pr=mdb.models[ModelName].rootAssembly.getMassProperties()
	mass=pr['mass']
	print (mass)

	print 'Made it in'
	odbName = ModelName+'.odb'
	print odbName
	odb = visualization.openOdb(odbName)
	lastFrame = odb.steps['Step-1'].frames[-1]
	print 'odb open'

	## The following is left for use in later probems/projects
	print 'Scanning the PART for maximum VM STRESS'
	elsetName='ALL_PART'
	elset = elemset = None
	region = "over the entire model"
	assembly = odb.rootAssembly

	lastFrame = odb.steps['Step-1'].frames[-1]
	print 'odb open'
 
	# Select Node(s) To Be Queried
	pTip = odb.rootAssembly.nodeSets['TIPNODE']
	# 604 ONLY Change for Optimization:
	# pTip = odb.rootAssembly.instances['PART-1-1'].nodeSets['TIPNODE']

	# Retrieve Y-displacements at Splines/Connectors
	print 'Retrieving ALL final displacements at ALL points'
	dispField = lastFrame.fieldOutputs['U']

	print 'Retrieving ALL displacements at TIPNODE'
	dFieldpTip = dispField.getSubset(region=pTip)
	
	print 'Retrieving only U2 at TIPNODE'
	# Note that U1=data[0], U2=data[1], U3=data[2]
	tipDisp = dFieldpTip.values[0].data[1]

	## The following is left for use in later probems/projects
	print 'Scanning the PART for maximum VM STRESS'
	elsetName='ALL_PART'
	elset = elemset = None
	region = "over the entire model"
	assembly = odb.rootAssembly
	
	if elsetName:
		try:
			elemset = assembly.elementSets[elsetName]
			region = " in the element set : " + elsetName;
		except KeyError:
			print 'An assembly level elset named %s does' \
				  'not exist in the output database %s' \
				  % (elsetName, odbName)			
			
	""" Initialize maximum values """
	maxMises = -0.1
	maxVMElem = 0
	maxStep = "_None_"
	maxFrame = -1
	Stress = 'S'
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
	print(type(maxMises))
	print(type(mass))
	print(type(volume))
	odb.close()

	Mass = volume*2720
	
	
	
	toc = time.time()
	elapsedTime = toc - tic
	print elapsedTime

	DataFile = open(filename, 'a')
	DataFile.write('%10f %10f %10f\n' %(maxMises, tipDisp, Mass))
	# DataFile.write('elapsed time = %10f seconds\n' % elapsedTime)
	DataFile.close()

print 'DONE!!'

