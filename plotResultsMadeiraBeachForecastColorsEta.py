##%% Plot TWL estimates on raw camera image
##% also plots measured runup values on camera image
##% plots cross-shore profile with wl values along it

### IMPORTS ###
import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import numpy as np


### FUNCTIONS ###
def coordSys_madbeach(E, N):
    '''
    Transform UTM to a local "madbeach" coordinate system

     Input:
       E, N (float) = UTM coordinates

     Output:
       X, Y (float) = local coordinates

     Usage: x, y= coordSys_madbeach(E,N)

     Written By:
       Jenna Brown, 2017, USGS
           based on Dave Thompson's madbeach_transformUTMToLocal.m
    '''

    #parameters from Madeira Beach calibration (intrinsics, extrinsics) - 01/28/2022
    #angle of rotation, degrees
    theta = 227.753111137952
    #origin Easting
    E0 = 323055.06675
    #origin Northing
    N0 = 3075922.21413775

    xyRotate(E, N, theta, xo=E0, yo=N0)

def xyRotate(x, y, theta, xo=0, yo=0):
    '''
    xyRotate  Rotate data.
    XR, YR = xyRotate(x, y, theta) rotates the coordinates x, y, and theta
    (cartesian decimal degrees).

    XR, YR = xyRotate(x, y, theta, xo, yo) rotates the coordinates around the
    origin xo, yo.

    Inputs:
        x (float) - x coordinate
        y (float) - y coordinate
        theta (float) - angle of rotation
        xo (float) - x origin
        yo (float) - y origin
    Outputs:
        XR (float) - rotated x coordinate
        YR (float) - rotated y coordinate
    '''
    #convert theta to radians
    theta = np.deg2rad(theta)
    
    #rotation matrix. Convert theta from degrees to radian for correct output
    A = np.array([[np.cos(theta), np.sin(theta)], [(-1)*np.sin(theta), np.cos(theta)]])

    columnE = x.transpose()
    print(columnE)
    print(columnE.shape)

    
    

### MAIN ###
#load image data
snapFile = '1665860400.Sat.Oct.15_19_00_00.GMT.2022.madbeach.c1.snap.jpg'
timexFile = '1665860400.Sat.Oct.15_19_00_00.GMT.2022.madbeach.c1.timex.jpg'
snap = plt.imread(snapFile)
timex = plt.imread(timexFile)
geom = scipy.io.loadmat('./matlabcode/geomFile_c1.mat')

#get datetime elements from filename
filenameElements = snapFile.split('.')
epoch = filenameElements[0]
dayOfWeek = filenameElements[1]
shortMonth = filenameElements[2]
dayHourSecond = filenameElements[3]
day = dayHourSecond.split('_')[0]
hour = dayHourSecond.split('_')[1]
second = dayHourSecond.split('_')[2]
year = filenameElements[5]

###display image
##plt.imshow(snap)
##plt.show()
##plt.imshow(timex)
##plt.show()

surveysAll = os.listdir('//gs/stpetersburgfl-g/NACCH/Imagery/madbeach/surveys/walking/')
i = 0
for elem in surveysAll:
    #want only most recent folder
    if elem.startswith('202'):
        if i == 0:
            mostRecent = elem
            i = i + 1
        else:
            if int(elem) > int(mostRecent):
                #2nd most recent in case most recent survey hasn't been processed yet
                secondMostRecent = mostRecent
                mostRecent = elem
                

mostRecentDir = os.listdir('//gs/stpetersburgfl-g/NACCH/Imagery/madbeach/surveys/walking/' + mostRecent)
#check that most recent survey folder has .xyz files
hasXYZ = False
for file in mostRecentDir:
    if file.endswith('.xyz'):
        hasXYZ = True
        break

#yes, most recent folder has .xyz files
if hasXYZ:
    mostRecentLine60 = '//gs/stpetersburgfl-g/NACCH/Imagery/madbeach/surveys/walking/' + mostRecent + '/line60.xyz'
else:
    mostRecentLine60 = '//gs/stpetersburgfl-g/NACCH/Imagery/madbeach/surveys/walking/' + secondMostRecent + '/line60.xyz'

E = []
N = []
with open(mostRecentLine60, 'r') as xyzFile:
    for line in xyzFile:
        xyz = line.split(' ')
        x = float(xyz[0])
        y = float(xyz[1])
        z = xyz[2]
        #z has newiine character at the end
        z = float(z.replace('\n', ''))
        E.append(x)
        N.append(y)

E = np.array(E)
N = np.array(N)
coordSys_madbeach(E, N)
