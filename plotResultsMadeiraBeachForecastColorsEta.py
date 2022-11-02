##%% Plot TWL estimates on raw camera image
##% also plots measured runup values on camera image
##% plots cross-shore profile with wl values along it

### IMPORTS ###
import os
import datetime
import time
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
from numpy import random

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

    X,Y = xyRotate(E, N, theta, xo=E0, yo=N0)

    return X,Y

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
    
    outX = x - xo
    outY = y - yo
    out = np.array([outX, outY])
    #transpose from shape 2,len(247) to shape len(x),247. e.g. 2,247 --> 247,2
    out = np.transpose(out)
    out = np.matmul(out, A)

    XR = out[0:len(x), 0]
    YR = out[0:len(y), 1]

    return XR,YR
    
    

### MAIN ###
#load image data
snapFile = '1663245000.Thu.Sep.15_12_30_00.GMT.2022.madbeach.c1.snap.jpg'
timexFile = '1663245000.Thu.Sep.15_12_30_00.GMT.2022.madbeach.c1.timex.jpg'
snap = plt.imread(snapFile)
timex = plt.imread(timexFile)
geom_c1 = scipy.io.loadmat('./matlabcode/geomFile_c1.mat')
geom_c2 = scipy.io.loadmat('./matlabcode/geomFile_c2.mat')

#get datetime elements from filename
filenameElements = snapFile.split('.')
epoch = filenameElements[0]
dayOfWeek = filenameElements[1]
shortMonth = filenameElements[2]
month = time.strptime(shortMonth, '%b').tm_mon
dayHourSecond = filenameElements[3]
day = dayHourSecond.split('_')[0]
hour = dayHourSecond.split('_')[1]
minute = dayHourSecond.split('_')[2]
second = dayHourSecond.split('_')[3]
year = filenameElements[5]
imgDatetimeStr = f'{year}-{month}-{day} {hour}:{minute}:{second}'
imgDatetime = datetime.datetime.strptime(imgDatetimeStr, "%Y-%m-%d %H:%M:%S")

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

#process line60
profileE = []
profileN = []
profileZ = []
with open(mostRecentLine60, 'r') as xyzFile:
    for line in xyzFile:
        xyz = line.split(' ')
        x = float(xyz[0])
        y = float(xyz[1])
        z = xyz[2]
        #z has newiine character at the end
        z = float(z.replace('\n', ''))
        profileE.append(x)
        profileN.append(y)
        profileZ.append(z)

profileE = np.array(profileE)
profileN = np.array(profileN)
profileZ = np.array(profileZ)
profileX, profileY = coordSys_madbeach(profileE, profileN)

#change NaN to 0s
for i, elem in enumerate(profileZ):
    if (elem == np.nan) or (elem == None):
        profileZ[i] = 0

uniqueZ = np.unique(profileZ)
a, b = np.histogram(profileZ, bins=uniqueZ)
mult = np.argwhere(a > 1)
#catch >2 identical values in profileZ
for i in range(0, len(mult)):
    #b is ascending array of unique values. mult is array of bins (arrays) in a with dupllicate values. indices in a match indices in b. So, index mult[i][0] will give index of a value in b that is not unique
    dupValue = b[mult[i][0]]
    dupIndices = np.argwhere(profileZ == dupValue)
    
    for ii in range(0, len(dupIndices)):
        #dupIndices[ii][0] is index in profileZ that contains duplicate values
        index = dupIndices[ii][0]
        #slightly alter value to prevent duplicates
        offset = random.rand() * 0.0000000001
        profileZ[index] = profileZ[index] + offset

uniqueZ = np.unique(profileZ)
a, b = np.histogram(profileZ, bins=uniqueZ)
mult = np.argwhere(a > 1)

lat = 27.796216949206798
lon = -82.796102542950635

all_twl = scipy.io.loadmat('./matlabcode/all_TWL_forecast.mat')
Rtime = all_twl['all_twl']['forecastTime']
Rrunup05 = all_twl['all_twl']['runup05']
Rrunup = all_twl['all_twl']['runup']
Rrunup95 = all_twl['all_twl']['runup95']
Rtwl05 = all_twl['all_twl']['twl05']
Rtwl = all_twl['all_twl']['twl']
Rtwl95 = all_twl['all_twl']['twl95']
Rslope = all_twl['all_twl']['slope']

#only save madbeach data
pindex = 1
#find twl data that matches image date. Because of scipy loadmat(), array is 3-layer deep array
for i in range(0, len(Rtime)):
    #convert string to datetime object
    RDatetime = datetime.datetime.strptime(Rtime[i][0][0], "%Y-%m-%d %H:%M:%S.0")

    if (RDatetime < imgDatetime + datetime.timedelta(hours=1)) and (RDatetime > imgDatetime - datetime.timedelta(hours=1)):
        print(RDatetime)

