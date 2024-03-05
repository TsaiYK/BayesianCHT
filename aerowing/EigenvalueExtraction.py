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

def findEigenValue(ModelName,StepName):
    odbName = ModelName+'.odb'
    odb = visualization.openOdb(odbName)

    lastFrameM = odb.steps[StepName].frames[-1]
    print lastFrameM.description
##The following string is from the description of the frame; contains the eigenvalue
    descString=lastFrameM.description
    print lastFrameM.mode
##Now we split the string at the = sign
    pattern2 = re.compile('\s*=\s*')
    print pattern2.split(descString)
    print pattern2.split(descString)[1]
##Convert the second string (index=1) to floating point number
    eigenVal1=float(pattern2.split(descString)[1])
##Test that it is a number and not a string
#    print eigenVal1/20.
    odb.close()

    return eigenVal1